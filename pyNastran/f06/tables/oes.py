#pylint: disable=C0301,C0111
from __future__ import print_function
from six.moves import range
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStress, RealRodStrain
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStress, RealBarStrain
#from pyNastran.op2.tables.oes_stressStrain.real.oes_beams   import RealBeamStress
#from pyNastran.op2.tables.oes_stressStrain.real.oes_shear   import RealShearStress
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStress, RealSolidStrain
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStress, RealPlateStrain
from pyNastran.op2.tables.oes_stressStrain.real.oes_compositePlates import RealCompositePlateStress, RealCompositePlateStrain
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import RealCelasStress, RealCelasStrain
#strain...


class OES(object):
    def _parse_line_blanks(self, sline, data_types, debug=False):
        pass

    def __init__(self):
        self.log = None
        self.iSubcases = []

    def _stress_in_crod_elements(self):
        self._stress_in_rod_elements('CROD', 1, 'crod_stress')

    def _stress_in_ctube_elements(self):
        self._stress_in_rod_elements('CTUBE', 3, 'ctube_stress')

    def _stress_in_conrod_elements(self):
        self._stress_in_rod_elements('CONROD', 10, 'crod_stress')

    def _strain_in_crod_elements(self):
        self._strain_in_rod_elements('CROD', 1, 'crod_strain')

    def _strain_in_ctube_elements(self):
        self._strain_in_rod_elements('CTUBE', 3, 'ctube_strain')

    def _strain_in_conrod_elements(self):
        self._strain_in_rod_elements('CONROD', 10, 'crod_stress')

    def _stress_in_rod_elements(self, element_name, element_type, result_type):
        """
        ::

                                       S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
          ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
            ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
               14    2.514247E+04              1.758725E+02                     15    2.443757E+04              2.924619E+01
        """
        (isubcase, transient, dt, data_code) = self._get_rod_header(element_name, element_type, False)
        data_code['table_name'] = 'OES1X'
        data = self._read_rod_stress()
        slot = getattr(self, result_type)
        if isubcase in slot:
            slot[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            slot[isubcase] = RealRodStress(data_code, is_sort1, isubcase, transient)
            slot[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strain_in_rod_elements(self, element_name, element_type, result_name):
        slot = getattr(self, result_name)
        (isubcase, transient, dt, data_code) = self._get_rod_header(element_name, element_type, False)
        data_code['table_name'] = 'OSTR1X'
        data = self._read_rod_stress()
        if isubcase in slot:
            slot[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            slot[isubcase] = RealRodStrain(data_code, is_sort1, isubcase, transient)
            slot[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)


    def _stress_in_celas1_elements(self):
        self._stress_in_celas_elements('CELAS2', 11, 'celas1_stress')

    def _stress_in_celas2_elements(self):
        self._stress_in_celas_elements('CELAS2', 12, 'celas2_stress')

    def _stress_in_celas3_elements(self):
        self._stress_in_celas_elements('CELAS3', 13, 'celas3_stress')

    def _stress_in_celas4_elements(self):
        self._stress_in_celas_elements('CELAS4', 14, 'celas4_stress')

    def _stress_in_celas_elements(self, element_name, element_type, result_name):
        slot = getattr(self, result_name)
        (isubcase, transient, dt, data_code) = self._get_spring_header(element_name, element_type, False)
        data_code['table_name'] = 'OSTR1X'
        data = self._read_spring_stress()
        if isubcase in slot:
            slot[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            slot[isubcase] = RealCelasStress(data_code, is_sort1, isubcase, transient)
            slot[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strain_in_celas1_elements(self):
        self._strain_in_celas_elements('CELAS2', 11, 'celas1_strain')

    def _strain_in_celas2_elements(self):
        self._strain_in_celas_elements('CELAS2', 12, 'celas2_strain')

    def _strain_in_celas3_elements(self):
        self._strain_in_celas_elements('CELAS3', 13, 'celas3_strain')

    def _strain_in_celas4_elements(self):
        self._strain_in_celas_elements('CELAS4', 14, 'celas4_strain')

    def _strain_in_celas_elements(self, element_name, element_type, result_name):
        slot = getattr(self, result_name)
        (isubcase, transient, dt, data_code) = self._get_spring_header(element_name, element_type, True)
        data_code['table_name'] = 'OSTR1X'
        data = self._read_spring_stress()
        if isubcase in slot:
            slot[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            slot[isubcase] = RealCelasStrain(data_code, is_sort1, isubcase, transient)
            slot[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _get_rod_header(self, element_name, element_type, is_strain):
        """
        * analysis_code = 1 (Statics)
        * device_code   = 1 (Print)
        * table_code    = 5 (Stress)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (Real)
        * s_code        = 0 (Stress)
        * num_wide      = 8 (???)
        """
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)

        (stress_bits, s_code) = make_stress_bits(is_strain=False, is_rod_or_solid=True)
        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 5, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
            'stress_bits': stress_bits, 'format_code': 1,
            'element_name': element_name, 'element_type': element_type, 'nonlinear_factor': dt,
            'lsdvmn' : 1,
            'dataNames':['lsdvmn']
        }

        return (isubcase, transient, dt, data_code)

    def _read_rod_stress(self):
        """
        ::

                                       S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
          ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
            ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
               14    2.514247E+04              1.758725E+02                     15    2.443757E+04              2.924619E+01
        """
        data = []
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            sline = [line[0:13], line[13:29], line[29:42], line[42:55], line[55:67], line[67:78], line[78:94], line[94:107], line[107:120], line[120:131]]
            if 'PAGE' in line:
                break
            #print sline
            data_types = [int, float, float, float, float, int, float, float, float, float]
            out = self._parse_line_blanks(sline, data_types)  # line 1
            #print out
            data.append(out[:5])
            if isinstance(out[5], int):
                data.append(out[5:])
            self.i += 1
        return data

    def _read_spring_stress(self):
        """
        ::

                                     S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 2 )
            ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN           ELEMENT         STRAIN
              ID.                              ID.                              ID.                              ID.
              20001      0.0                   20002      0.0                   20003      0.0                   20004      0.0
              20005      0.0                   20006      0.0
         """
        data = []
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            sline = line.strip().split()
            #sline = [line[0:13], line[13:29], line[29:42], line[42:55], line[55:67], line[67:78], line[78:94], line[94:107], line[107:120], line[120:131]]
            if 'PAGE' in line:
                break
            #print sline
            n = len(sline) // 2
            assert len(sline) % 2 == 0, sline

            data_types = [int, float] * n
            out = self._parse_line_blanks(sline, data_types)  # line 1
            #print out

            while out:
                strain = out.pop()
                eid = out.pop()
                data.append([eid, strain])
            self.i += 1

        return data

    def _stress_in_cbar_elements(self):
        """
        ::

                                         S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )
          ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
            ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
               12    0.0            0.0            0.0            0.0            1.020730E+04   1.020730E+04   1.020730E+04
                     0.0            0.0            0.0            0.0                           1.020730E+04   1.020730E+04

        * analysis_code = 1 (Statics)
        * device_code   = 1 (Print)
        * table_code    = 5 (Stress)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (Real)
        * s_code        = 0 (Stress)
        * num_wide      = 8 (???)
        """
        (isubcase, transient, dt, data_code) = self._get_bar_header(False)

        data = self._read_bar_stress()
        data_code['table_name'] = 'OES1X'
        if isubcase in self.cbar_stress:
            self.cbar_stress[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.cbar_stress[isubcase] = RealBarStress(data_code, is_sort1, isubcase, dt)
            self.cbar_stress[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strain_in_cbar_elements(self):
        (isubcase, transient, dt, data_code) = self._get_bar_header(False)
        data_code['table_name'] = 'OSTR1X'

        data = self._read_bar_stress()
        if isubcase in self.cbar_strain:
            self.cbar_strain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.cbar_strain[isubcase] = RealBarStrain(data_code, is_sort1, isubcase, dt)
            self.cbar_strain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _get_bar_header(self, is_strain):
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        (stress_bits, s_code) = make_stress_bits(is_strain=is_strain, is_rod_or_solid=True)
        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 5, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
            'stress_bits': stress_bits, 'format_code': 1,
            'element_name': 'CBAR', 'element_type': 34,
            'nonlinear_factor': dt,
            'lsdvmn' : 1,
            'dataNames':['lsdvmn']
        }

        return (isubcase, transient, dt, data_code)

    def _get_spring_header(self, element_name, element_type, is_strain):
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        (stress_bits, s_code) = make_stress_bits(is_strain=is_strain, is_rod_or_solid=False)
        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 5, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 2, 's_code': s_code,
            'thermal': 0,
            'stress_bits': stress_bits, 'format_code': 1,
            'element_name': element_name, 'element_type': element_type,
            'nonlinear_factor': dt,
            'dataNames':['lsdvmn']
        }
        return (isubcase, transient, dt, data_code)

    def _read_bar_stress(self):
        """
        ::

          ELEMENT        SA1            SA2            SA3            SA4           AXIAL          SA-MAX         SA-MIN     M.S.-T
            ID.          SB1            SB2            SB3            SB4           STRESS         SB-MAX         SB-MIN     M.S.-C
               12    0.0            0.0            0.0            0.0            1.020730E+04   1.020730E+04   1.020730E+04
                     0.0            0.0            0.0            0.0                           1.020730E+04   1.020730E+04
        """
        data = []
        while 1:
            # line 1
            line = self.infile.readline()[1:].rstrip('\r\n ')
            if 'PAGE' in line:
                break
            sline = [line[0:11], line[11:26], line[26:41], line[41:56], line[56:69], line[69:86], line[86:101], line[101:115], line[115:131]]
            data_types = [int, float, float, float, float, float, float, float, float]
            out = self._parse_line_blanks(sline, data_types)  # line 1
            out = ['CBAR'] + out
            #data.append(sline)

            # line 2
            line = self.infile.readline()[1:].rstrip('\r\n ')
            sline = [line[11:26], line[26:41], line[41:56], line[56:69], line[86:101], line[101:115], line[115:131]]
            data_types = [float, float, float, float, float, float, float]
            out += self._parse_line_blanks(sline, data_types)  # line 2
            #print "*",out
            data.append(out)
            self.i += 2

        return data

    #==========================================================================
    # COMPOSITE
    def _stress_in_composite_ctria3_elements(self):
        element_name = 'CTRIA3'
        element_type = 97
        is_strain = False
        self._composites_helper(element_name, element_type, is_strain, self.ctria3_composite_stress)

    def _strain_in_composite_ctria3_elements(self):
        element_name = 'CTRIA3'
        element_type = 97
        is_strain = True
        self._composites_helper(element_name, element_type, is_strain, self.ctria3_composite_strain)

    def _stress_in_composite_cquad4_elements(self):
        """
        ::

                         S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
          ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
            ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
              181    1   3.18013E+04  5.33449E+05  1.01480E+03   -7.06668E+01  1.90232E+04   89.88  5.33451E+05  3.17993E+04  2.50826E+05
              181    2   1.41820E+05  1.40805E+05  1.25412E+05   -1.06000E+02  2.85348E+04   44.88  2.66726E+05  1.58996E+04  1.25413E+05

        element_type = 33 b/c not bilinear
        """
        element_name = 'CQUAD4'
        element_type = 95
        is_strain = False
        self._composites_helper(element_name, element_type, is_strain, self.cquad4_composite_stress)

    def _strain_in_composite_cquad4_elements(self):
        """
        ::

                         S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )
          ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS    INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)      MAX
            ID      ID    NORMAL-1     NORMAL-2     SHEAR-12     SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE    MAJOR        MINOR        SHEAR
              181    1   3.18013E+04  5.33449E+05  1.01480E+03   -7.06668E+01  1.90232E+04   89.88  5.33451E+05  3.17993E+04  2.50826E+05
              181    2   1.41820E+05  1.40805E+05  1.25412E+05   -1.06000E+02  2.85348E+04   44.88  2.66726E+05  1.58996E+04  1.25413E+05

        element_type = 33 b/c not bilinear
        """
        element_name = 'CQUAD4'
        element_type = 95
        is_strain = True
        self._composites_helper(element_name, element_type, is_strain, self.cquad4_composite_strain)

    def _composites_helper(self, element_name, element_type, is_strain, slot):
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        data_types = [int, int, float, float, float, float,
                      float, float, float, float, float]
        data = self._read_f06_table(data_types)

        is_max_shear = False  # Von Mises/Max Shear
        sHeaders = headers.rstrip()
        if 'SHEAR' in sHeaders[-5:]:  # last 5 letters of the line to avoid 'SHEAR YZ-MAT'
            is_max_shear = True
        else:
            raise RuntimeError(sHeader)

        is_fiber_distance = True
        (stress_bits, s_code) = make_stress_bits(is_fiber_distance=is_fiber_distance, is_max_shear=is_max_shear, is_strain=is_strain)
        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 5, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
            'stress_bits': stress_bits, 'format_code': 1,
            'element_name': element_name, 'element_type': element_type,
            'table_name': 'OES1X', 'nonlinear_factor': dt,
            'dataNames':['lsdvmn'],
            'lsdvmn': 1,
        }

        is_sort1 = True
        if is_strain:
            class_obj = RealCompositePlateStrain
        else:
            class_obj = RealCompositePlateStress

        if isubcase in slot:
            slot[isubcase].add_f06_data(data, transient, element_name)
        else:
            assert 'nonlinear_factor' in data_code
            slot[isubcase] = class_obj(data_code, is_sort1, isubcase, transient)
            slot[isubcase].add_f06_data(data, transient, element_name)
        self.iSubcases.append(isubcase)
    #==========================================================================

    def _stress_in_ctria3_elements(self):
        """
        ::

                                   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )
          ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)
            ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
                8   -1.250000E-01     -1.303003E+02   1.042750E+04  -1.456123E+02   -89.2100    1.042951E+04   -1.323082E+02   1.049629E+04
                     1.250000E-01     -5.049646E+02   1.005266E+04  -2.132942E+02   -88.8431    1.005697E+04   -5.092719E+02   1.032103E+04

        * analysis_code = 1 (Statics)
        * device_code   = 1 (Print)
        * table_code    = 5 (Stress)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (Real)
        * s_code        = 0 (Stress)
        * num_wide      = 8 (???)
        """
        (isubcase, transient, data_code) = self._get_tri_header(False)
        data_code['table_name'] = 'OES1X'
        data = self._read_tri_stress(['CTRIA3'])
        if isubcase in self.ctria3_stress:
            self.ctria3_stress[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.ctria3_stress[isubcase] = RealPlateStress(data_code, is_sort1, isubcase, transient)
            self.ctria3_stress[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strain_in_ctria3_elements(self):
        (isubcase, transient, data_code) = self._get_tri_header(True)
        data_code['table_name'] = 'OST1X'
        data = self._read_tri_stress(['CTRIA3'])
        if isubcase in self.ctria3_strain:
            self.ctria3_strain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            assert 'nonlinear_factor' in data_code
            self.ctria3_strain[isubcase] = RealPlateStrain(data_code, is_sort1, isubcase, transient)
            self.ctria3_strain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _get_tri_header(self, is_strain):
        """
        * analysis_code = 1 (Statics)
        * device_code   = 1 (Print)
        * table_code    = 5 (Stress)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (Real)
        * s_code        = 0 (Stress)
        * num_wide      = 8 (???)
        """
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        is_fiber_distance = None
        is_max_shear = None  # Von Mises/Max Shear
        if 'DISTANCE' in headers:
            is_fiber_distance = True
        elif 'CURVATURE' in headers:
            is_fiber_distance = False
        else:
            raise RuntimeError(headers)

        if 'VON MISES' in headers:
            is_max_shear = False
        elif 'SHEAR' in headers:
            is_max_shear = True
        else:
            raise RuntimeError(headers)

        (stress_bits, s_code) = make_stress_bits(is_fiber_distance, is_max_shear, is_strain=is_strain)
        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 5, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
            'stress_bits': stress_bits, 'format_code': 1,
            'element_name': 'CTRIA3', 'element_type': 74,
            'nonlinear_factor': dt,
            'dataNames':['lsdvmn'],
            'lsdvmn': 1,
            }
        if transient is not None:
            data_code['name'] = transient[0]
        return (isubcase, transient, data_code)

    def _read_tri_stress(self, eType):
        """
        ::

                  ELEMENT      FIBER               STRESSES IN ELEMENT COORD SYSTEM             PRINCIPAL STRESSES (ZERO SHEAR)
                    ID.       DISTANCE           NORMAL-X       NORMAL-Y      SHEAR-XY       ANGLE         MAJOR           MINOR        VON MISES
                        8   -1.250000E-01     -1.303003E+02   1.042750E+04  -1.456123E+02   -89.2100    1.042951E+04   -1.323082E+02   1.049629E+04
                             1.250000E-01     -5.049646E+02   1.005266E+04  -2.132942E+02   -88.8431    1.005697E+04   -5.092719E+02   1.032103E+04
        """
        data = []
        while 1:
            line = self.infile.readline()[1:].strip().split()
            if 'PAGE' in line or '***' in line:
                break
            #print line
            data_types = [int, float, float, float, float,
                          float, float, float, float]
            sline = self.parseLine(line, data_types)  # line 1
            #print 'sline',sline
            sline = eType + sline
            data.append(sline)
            line = self.infile.readline()[1:].strip().split()
            #print '***',line
            data_types = [float, float, float, float,
                          float, float, float, float]
            sline += self.parseLine(line, data_types)  # line 2
            data.append(sline)
            self.i += 2
        return data

    #==========================================================================
    # CQUAD4
    def _stress_in_cquad4_elements(self):
        (isubcase, transient, data_code) = self._get_quad_header(2, 'CQUAD4', 33, is_strain=False)
        data_code['table_name'] = 'OES1X'
        data = self._read_tri_stress(['CQUAD4'])
        if isubcase in self.cquad4_stress:
            self.cquad4_stress[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            assert 'nonlinear_factor' in data_code
            self.cquad4_stress[isubcase] = RealPlateStress(data_code, is_sort1, isubcase, transient)
            self.cquad4_stress[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strain_in_cquad4_elements(self):
        (isubcase, transient, data_code) = self._get_quad_header(2, 'CQUAD4', 33, is_strain=True)
        data_code['table_name'] = 'OSTR1X'
        data = self._read_tri_stress(['CQUAD4'])
        if isubcase in self.cquad4_strain:
            self.cquad4_strain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.cquad4_strain[isubcase] = RealPlateStrain(data_code, is_sort1, isubcase, transient)
            self.cquad4_strain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strain_in_cquad4_bilinear_elements(self):
        self._stress_strain_cquad4_bilinear_helper('CQUAD4', 144, is_strain=True)

    def _stress_in_cquad4_bilinear_elements(self):
        self._stress_strain_cquad4_bilinear_helper('CQUAD4', 144, is_strain=False)

    def _stress_strain_cquad4_bilinear_helper(self, element_type, element_num, is_strain=None):
        """
        ::

                               S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN

          ELEMENT              FIBER            STRESSES IN ELEMENT COORD SYSTEM         PRINCIPAL STRESSES (ZERO SHEAR)
            ID      GRID-ID   DISTANCE        NORMAL-X      NORMAL-Y      SHEAR-XY      ANGLE        MAJOR         MINOR       VON MISES
                6    CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  8.247786E+03
                             1.250000E-01   5.406062E+02  1.201854E+04 -4.174177E+01   -89.7916   1.201869E+04  5.404544E+02  1.175778E+04

                         4  -1.250000E-01  -8.871141E+02  7.576036E+03 -1.550089E+02   -88.9511   7.578874E+03 -8.899523E+02  8.060780E+03
                             1.250000E-01  -8.924081E+01  1.187899E+04 -4.174177E+01   -89.8002   1.187913E+04 -8.938638E+01  1.192408E+04
        """
        (isubcase, transient, data_code, is_max_shear) = self._get_quad_header(3, element_type, element_num, is_strain)
        #print(self.getQuadHeader(3, False, 144))
        #print("data_code =", data_code)

        data_code['table_name'] = 'OES1X'
        data = self._read_quad_bilinear()

        if is_strain:
            slot = self.cquad4_strain
            class_obj = RealPlateStrain
        else:
            slot = self.cquad4_stress
            class_obj = RealPlateStress


        if isubcase in slot:
            slot[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            assert 'nonlinear_factor' in data_code
            result = class_obj(data_code, is_sort1, isubcase, transient)
            result.add_f06_data(data, transient)
            slot[isubcase] = result
            is_von_mises = not(is_max_shear)
            assert result.is_max_shear() == is_max_shear
            assert result.is_von_mises() == is_von_mises
        self.iSubcases.append(isubcase)

    def _get_quad_header(self, nHeaderLines, element_name, element_type, is_strain):
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(nHeaderLines)
        #print "headers = %s" %(headers)

        is_fiber_distance = None
        is_max_shear = None  # Von Mises/Max Shear
        if 'DISTANCE' in headers:
            is_fiber_distance = True
        elif 'CURVATURE' in headers:
            is_fiber_distance = False
        else:
            raise RuntimeError(headers)

        if 'SHEAR' in headers:
            is_max_shear = True
        elif 'VON MISES' in headers:
            is_max_shear = False
        else:
            raise RuntimeError(headers)

        (stress_bits, s_code) = make_stress_bits(is_fiber_distance, is_max_shear, is_strain)
        data_code = {
            'analysis_code': analysis_code,
            'device_code': 1, 'table_code': 5, 'sort_code': 0,
            'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
            'stress_bits': stress_bits, 'format_code': 1,
            'element_name': element_name, 'element_type': element_type,
            'nonlinear_factor': dt,
            'dataNames':['lsdvmn'],
            'lsdvmn': 1,
        }
        if transient is not None:
            data_code['name'] = transient[0]
        return (isubcase, transient, data_code, is_max_shear)

    def _read_quad_bilinear(self):
        data = []
        while 1:
            if 1:  # CEN/4
                line = self.infile.readline()[1:].strip().split()
                if 'PAGE' in line:
                    return data
                sline = self.parseLine(line, [int, str, float, float, float, float, float, float, float, float])  # line 1
                if sline is None:
                    break
                sline = ['CQUAD4'] + sline
                #data.append(sline)
                line = self.infile.readline()[1:].strip().split()
                data_types = [float, float, float, float,
                              float, float, float, float]
                sline += self.parseLine(line, data_types)  # line 2
                data.append(sline)
                line = self.infile.readline()  # blank line
                self.i += 3

            for i in range(4):
                line = self.infile.readline()[1:].strip().split()
                data_types = [int, float, float, float, float,
                              float, float, float, float]
                sline = self.parseLine(line, data_types)  # line 1
                #data.append(sline)
                line = self.infile.readline()[1:].strip().split()
                data_types = [float, float, float, float,
                              float, float, float, float]
                sline += self.parseLine(line, data_types)  # line 2
                data.append(sline)
                line = self.infile.readline()  # blank line
                self.i += 3

        return data

    #==========================================================================

    def _stress_in_chexa_elements(self):
        return self._read_solid_stress('CHEXA', 67, 8, self.chexa_stress)

    def _stress_in_cpenta_elements(self):
        return self._read_solid_stress('CPENTA', 68, 6, self.cpenta_stress)

    def _stress_in_ctetra_elements(self):
        return self._read_solid_stress('CTETRA', 39, 4, self.ctetra_stress)

    def _strain_in_chexa_elements(self):
        return self._read_solid_strain('CHEXA', 67, 8, self.chexa_strain)

    def _strain_in_cpenta_elements(self):
        return self._read_solid_strain('CPENTA', 68, 6, self.cpenta_strain)

    def _strain_in_ctetra_elements(self):
        return self._read_solid_strain('CTETRA', 39, 4, self.ctetra_strain)

    def _read_solid_stress(self, element_name, element_type, nnodes, slot):
        (isubcase, transient, data_code) = self._get_solid_header(element_name, element_type, nnodes, is_strain=False)
        data_code['table_name'] = 'OES1X'

        data = self._read_3D_stress(element_name, nnodes, is_strain=False)
        if isubcase in slot:
            slot[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            slot[isubcase] = RealSolidStress(data_code, is_sort1, isubcase, transient)
            slot[isubcase].add_f06_data(data, transient)
            assert slot[isubcase].is_stress() == True
            assert slot[isubcase].is_strain() == False
        self.iSubcases.append(isubcase)

    def _read_solid_strain(self, element_name, element_type, nnodes, slot):
        (isubcase, transient, data_code) = self._get_solid_header(element_name, element_type, nnodes, is_strain=True)
        data_code['table_name'] = 'OSTR1X'

        data = self._read_3D_stress(element_name, nnodes, is_strain=True)
        if isubcase in slot:
            slot[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            slot[isubcase] = RealSolidStrain(data_code, is_sort1, isubcase, transient)
            slot[isubcase].add_f06_data(data, transient)
            assert slot[isubcase].is_stress() == False
            assert slot[isubcase].is_strain() == True
        self.iSubcases.append(isubcase)

    def _get_solid_header(self, element_name, element_type, n, is_strain=True):
        """
        * analysis_code = 1 (Statics)
        * device_code   = 1 (Print)
        * table_code    = 5 (Stress/Strain)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (Real)
        * s_code        = 0 (Stress/Strain)
        * num_wide      = 8 (???)
        """
        (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        headers = self.skip(2)
        #print("headers = %s" % (headers))

        is_max_shear = True
        if 'VON MISES' in headers:
            is_max_shear = False

        (stress_bits, s_code) = make_stress_bits(is_max_shear=False, is_strain=is_strain, is_rod_or_solid=True)
        data_code = {
            'analysis_code': 1, 'device_code': 1,
            'table_code': 5, 'sort_code': 0, 'sort_bits': [0, 0, 0],
            'num_wide': 8,
            'element_name': element_name, 'element_type': element_type,
            'format_code': 1,
            's_code': s_code, 'stress_bits': stress_bits,
            'nonlinear_factor': dt,
            'lsdvmn' : 1,
            'dataNames':['lsdvmn']
        }
        if transient is not None:
            data_code['name'] = transient[0]
        return (isubcase, transient, data_code)

    def _read_3D_stress(self, eType, nnodes, is_strain=False):
        data = []
        eType2 = eType + str(nnodes)
        while 1:
            line = self.infile.readline().rstrip('\n\r')  # [1:]
                    #              CENTER         X          #          XY             #        A         #
            sline = [line[1:17], line[17:23], line[23:28], line[28:43], line[43:47], line[47:63], line[63:66], line[66:80], line[80:83], line[83:88], line[88:93], line[93:98], line[99:113], line[113:130]]
            sline = [s.strip() for s in sline]
            if 'PAGE' in line:
                break
            elif '' is not sline[0]:
                sline = [eType2] + sline
            data.append(sline)

        return data

def make_stress_bits(is_fiber_distance=False, is_max_shear=True, is_strain=True, is_rod_or_solid=False):
    """
    Therefore, stress_code can be one of the following values:
    +------+---------+----------------------------------------------+
    |Value | On bits | Description                                  |
    +------+---------+----------------------------------------------+
    |  0   | 0 0 0 0 | Stress maximum shear or octahedral           |
    |  1   | 0 0 0 1 | Stress von Mises                             |
    |  10  | 1 0 1 0 | Strain Curvature maximum shear or octahedral |
    |  11  | 1 0 1 1 | Strain Curvature von Mises                   |
    |  14  | 1 1 1 0 | Strain Fibre maimum shear or octahedral      |
    |  15  | 1 1 1 1 | Strain Fibre von Mises                       |
    +------+---------+----------------------------------------------+
    """
    #print "is_max_shear=%s is_fiber_distance=%s" %(is_max_shear, is_fiber_distance)

   #code = (isVonMises, isFiberCurvatur, isStress, isNotRod)
    von_mises_code = 0 if is_max_shear else 1
    strain_code = 1 if is_strain else 0
    fiber_code = 1 if is_fiber_distance else 0
    mat_coord = 0
    if is_rod_or_solid:
        pass
    else:
        stress_bits = [mat_coord, strain_code, fiber_code, strain_code, von_mises_code]
        #fiber_code = 0
        #von_mises_code = 0
        #stress_bits = [1, 0, 0, 0, 0]
        #stress_bits.reverse()
        #print('bits=%s' % stress_bits)
        s_code = 0
        for i, codei in enumerate(reversed(stress_bits)):
            s_code += codei * 2**i
        #print('s_code=%s' % (s_code))
        return (stress_bits, s_code)

    # True, False, True, False
    #[zero, one, two, three, is_von_mises]
    code = (is_max_shear, is_fiber_distance, is_strain, is_rod_or_solid)
    mapper = {
        # element coordinate system (no material support)
        (True,  False, False, True) : ([0, 0, 0, 0, 0], 0),  # 0,  rod/csolid
        (False, False, False, True) : ([0, 0, 0, 0, 1], 1),  # 1,  rod/csolid strain
        (False, False, True,  True) : ([0, 1, 0, 1, 1], 1),   # ???, rod/csolid strain

        #(True,  False, True, False) : ([0, 1, 0, 1, 0], 10),  # 10
        #(False, False, True, False) : ([0, 1, 0, 1, 1], 11),  # 11

        #(True,  True, False, False) : ([0, 1, 1, 1, 0], 14),  # 14 - max shear, fiber_dist, stress, quad/tri
        #(True,  True, True,  False) : ([0, 1, 1, 1, 0], 14),  # 14 - max shear, fiber_dist, strain, quad/tri
        #(False, True, True,  False) : ([0, 1, 1, 1, 1], 15),  # 15 - von mises, fiber_dist, strain, quad/tri

        #(True,  False, False, False) : ([0, 0, 0, 0, 0], 0),  # 0,  composite
        #(False, True,  False, False) : ([0, 0, 0, 0, 1], 0),  # cquad4 bilinear ??? why do i need this...
    }
    try:
        (stress_bits, s_code) = mapper[code]
    except KeyError:
        msg = 'is_max_shear=%s is_fiber_distance=%s is_strain=%s is_rod_or_solid=%s' % (is_max_shear, is_fiber_distance, is_strain, is_rod_or_solid)
        raise RuntimeError(msg)

    #if is_max_shear==False:
    #    stress_bits[4] = 1 # Von Mises
    #if is_strain:
    #    #stress_bits[1] = stress_bits[3] = 1 # Strain
    #    stress_bits[1] = stress_bits[3] = 1 # Strain
    #if is_fiber_distance:
    #    stress_bits[2] = 1 # FiberDistance
    #print stress_bits
    #s_code = 0
    #for i,bit in enumerate(stress_bits):
    #    s_code += bit*2**i

    #def is_max_shear(self):
        #if self.stress_bits[4] == 0:
            #return True
        #return False

    #def is_curvature(self):
        #if self.s_code in [0, 1, 14, 15, 16, 17, 27, 30, 31]:  # fiber distance
            #return False
        #elif self.s_code in [10, 11, 26, ]:  # fiber curvature
            #return True
        #raise NotImplementedError('add s_code=%s' % self.s_code)

    #def isCurvatureOld(self):
        #if self.stress_bits[2] == 0:
            #return True
        #return False
    return (stress_bits, s_code)
