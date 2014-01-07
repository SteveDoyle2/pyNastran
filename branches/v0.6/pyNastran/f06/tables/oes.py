from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RodStressObject, RodStrainObject
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import BarStressObject, BarStrainObject
#from pyNastran.op2.tables.oes_stressStrain.real.oes_beams   import beamStressObject
#from pyNastran.op2.tables.oes_stressStrain.real.oes_shear   import shearStressObject
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import SolidStressObject, SolidStrainObject
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import PlateStressObject, PlateStrainObject
from pyNastran.op2.tables.oes_stressStrain.real.oes_compositePlates import CompositePlateStressObject, CompositePlateStrainObject
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import CelasStressObject, CelasStrainObject
#strain...


class OES(object):
    def __init__(self):
        self.rodStress = {}  # CROD, CONROD, CTUBE
        self.rodStrain = {}

        self.barStress = {}  # CBAR
        self.barStrain = {}

        self.plateStress = {}  # isotropic CTRIA3/CQUAD4
        self.plateStrain = {}

        self.solidStress = {}  # CTETRA/CPENTA/CHEXA
        self.solidStrain = {}

        self.compositePlateStress = {}  # composite CTRIA3/CQUAD4
        self.compositePlateStrain = {}

        #-------------
        # not supported
        self.celasStress = {}  # CELASi
        self.celasStrain = {}
        self.beamStress = {}  # CBEAM
        self.beamStrain = {}
        self.shearStress = {}  # CSHEAR
        self.shearStrain = {}
        self.nonlinearRodStress = {}  # CROD, CONROD, CTUBE
        self.nonlinearRodStrain = {}
        self.nonlinearPlateStress = {}  # CTRIA3, CTRIA6, CQUAD4, CQUAD8
        self.nonlinearPlateStrain = {}
        self.ctriaxStress = {}  # CTRIAX6
        self.ctriaxStrain = {}
        self.hyperelasticPlateStress = {}  # CTRIA3, CTRIA6, CQUAD4, CQUAD8
        self.hyperelasticPlateStrain = {}

    def _stresses_in_crod_elements(self):
        """
        ::

                                       S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )
          ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY       ELEMENT       AXIAL       SAFETY      TORSIONAL     SAFETY
            ID.        STRESS       MARGIN        STRESS      MARGIN         ID.        STRESS       MARGIN        STRESS      MARGIN
               14    2.514247E+04              1.758725E+02                     15    2.443757E+04              2.924619E+01
        """
        (isubcase, transient, dt, data_code) = self.getRodHeader(False)
        data_code['table_name'] = 'OES1X'
        data = self.readRodStress()
        if isubcase in self.rodStress:
            self.rodStress[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.rodStress[isubcase] = RodStressObject(data_code, is_sort1, isubcase, transient)
            self.rodStress[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strains_in_crod_elements(self):
        (isubcase, transient, dt, data_code) = self.getRodHeader(False)
        data_code['table_name'] = 'OSTR1X'
        data = self.readRodStress()
        if isubcase in self.rodStrain:
            self.rodStrain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.rodStrain[isubcase] = RodStrainObject(data_code, is_sort1, isubcase, transient)
            self.rodStrain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _stresses_in_celas2_elements(self):
        (isubcase, transient, dt, data_code) = self.getSpringHeader(False)
        data_code['table_name'] = 'OSTR1X'
        data = self.readSpringStress()
        if isubcase in self.celasStrain:
            self.celasStrain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.celasStrain[isubcase] = CelasStressObject(data_code, is_sort1, isubcase, transient)
            self.celasStrain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strains_in_celas2_elements(self):
        (isubcase, transient, dt, data_code) = self.getSpringHeader(False)
        data_code['table_name'] = 'OSTR1X'
        data = self.readSpringStress()
        if isubcase in self.celasStrain:
            self.celasStrain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.celasStrain[isubcase] = CelasStrainObject(data_code, is_sort1, isubcase, transient)
            self.celasStrain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def getRodHeader(self, isStrain):
        """
        * analysis_code = 1 (Statics)
        * device_code   = 1 (Print)
        * table_code    = 5 (Stress)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (Real)
        * s_code        = 0 (Stress)
        * num_wide      = 8 (???)
        """
        (subcaseName, isubcase, transient, dt, analysis_code, is_sort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        (stress_bits, s_code) = self.make_stress_bits(
            isStrain=False, isRodOrSolid=True)
        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 5, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
                    'stress_bits': stress_bits, 'format_code': 1,
                    'element_name': 'ROD', 'element_type': 1, 'nonlinear_factor': dt,
                    'dataNames':['lsdvmn']}

        return (isubcase, transient, dt, data_code)

    def readRodStress(self):
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
            dataTypes = [int, float, float, float, float, int, float, float, float, float]
            out = self.parseLineBlanks(sline, dataTypes)  # line 1
            #print out
            data.append(out[:5])
            if isinstance(out[5], int):
                data.append(out[5:])
            self.i += 1

        return data

    def readSpringStress(self):
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

            dataTypes = [int, float] * n
            out = self.parseLineBlanks(sline, dataTypes)  # line 1
            #print out

            while out:
                strain = out.pop()
                eid = out.pop()
                data.append([eid, strain])
            self.i += 1

        return data

    def _stresses_in_cbar_elements(self):
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
        (isubcase, transient, dt, data_code) = self.getBarHeader(False)

        data = self.readBarStress()
        data_code['table_name'] = 'OES1X'
        if isubcase in self.barStress:
            self.barStress[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.barStress[isubcase] = BarStressObject(data_code, is_sort1, isubcase, dt)
            self.barStress[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strains_in_cbar_elements(self):
        (isubcase, transient, dt, data_code) = self.getBarHeader(False)
        data_code['table_name'] = 'OSTR1X'

        data = self.readBarStress()
        if isubcase in self.barStrain:
            self.barStrain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.barStrain[isubcase] = BarStrainObject(data_code, is_sort1, isubcase, dt)
            self.barStrain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def getBarHeader(self, isStrain):
        (subcaseName, isubcase, transient, dt, analysis_code,
            is_sort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        (stress_bits, s_code) = self.make_stress_bits(isStrain=isStrain, isRodOrSolid=True)
        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 5, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
                    'stress_bits': stress_bits, 'format_code': 1,
                    'element_name': 'CBAR', 'element_type': 34,
                    'nonlinear_factor': dt,
                    'dataNames':['lsdvmn']}

        return (isubcase, transient, dt, data_code)

    def getSpringHeader(self, isStrain):
        (subcaseName, isubcase, transient, dt, analysis_code, is_sort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        (stress_bits, s_code) = self.make_stress_bits(isStrain=isStrain, isRodOrSolid=False)
        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 5, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 2, 's_code': s_code,
                    'thermal': 0,
                    'stress_bits': stress_bits, 'format_code': 1,
                    'element_name': 'ELAS2', 'element_type': 12,
                    'nonlinear_factor': dt,
                    'dataNames':['lsdvmn']}
        return (isubcase, transient, dt, data_code)

    def readBarStress(self):
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
            sline = [line[0:11], line[11:26], line[26:41], line[41:56], line[56:69], line[69:86], line[86:101], line[101:116], line[116:131]]
            dataTypes = [int, float, float, float, float, float, float, float, float]
            out = self.parseLineBlanks(sline, dataTypes)  # line 1
            out = ['CBAR'] + out
            #data.append(sline)

            # line 2
            line = self.infile.readline()[1:].rstrip('\r\n ')
            sline = [line[11:26], line[26:41], line[41:56], line[56:69], line[86:101], line[101:116], line[116:131]]
            dataTypes = [float, float, float, float, float, float, float]
            out += self.parseLineBlanks(sline, dataTypes)  # line 2
            #print "*",out
            data.append(out)
            self.i += 2

        return data

    #==========================================================================
    # COMPOSITE
    def _stresses_in_composite_ctria3_elements(self):
        element_name = 'CTRIA3'
        element_type = 97
        is_strain = False
        self._composites_helper(element_name, element_type, is_strain)

    def _strains_in_composite_ctria3_elements(self):
        element_name = 'CTRIA3'
        element_type = 97
        is_strain = True
        self._composites_helper(element_name, element_type, is_strain)

    def _stresses_in_composite_cquad4_elements(self):
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
        self._composites_helper(element_name, element_type, is_strain)

    def _strains_in_composite_cquad4_elements(self):
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
        self._composites_helper(element_name, element_type, is_strain)

    def _composites_helper(self, element_name, element_type, is_strain):
        (subcaseName, isubcase, transient, dt, analysis_code, is_sort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)
        dataTypes = [int, int, float, float, float, float,
                     float, float, float, float, float]
        data = self.readTable(dataTypes)

        isMaxShear = False  # Von Mises/Max Shear
        sHeaders = headers.rstrip()
        if 'SHEAR' in sHeaders[-5:]:  # last 5 letters of the line to avoid 'SHEAR YZ-MAT'
            isMaxShear = True
        else:
            raise RuntimeError(sHeader)

        isFiberDistance = True
        (stress_bits, s_code) = self.make_stress_bits(isFiberDistance=isFiberDistance, isMaxShear=isMaxShear, isStrain=is_strain)
        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 5, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
                    'stress_bits': stress_bits, 'format_code': 1,
                    'element_name': element_name, 'element_type': element_type,
                    'table_name': 'OES1X', 'nonlinear_factor': dt,
                    'dataNames':['lsdvmn']}

        if is_strain:
            dictA = self.compositePlateStrain
            class_obj = CompositePlateStrainObject
        else:
            dictA = self.compositePlateStress
            class_obj = CompositePlateStressObject

        if isubcase in dictA:
            dictA[isubcase].add_f06_data(data, transient, element_name)
        else:
            assert 'nonlinear_factor' in data_code
            dictA[isubcase] = class_obj(data_code, isubcase, transient)
            dictA[isubcase].add_f06_data(data, transient, element_name)
        self.iSubcases.append(isubcase)
    #==========================================================================

    def _stresses_in_ctria3_elements(self):
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
        (isubcase, transient, data_code) = self.getTriHeader(False)
        data_code['table_name'] = 'OES1X'
        data = self.readTriStress(['CTRIA3'])
        if isubcase in self.plateStress:
            self.plateStress[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.plateStress[isubcase] = PlateStressObject(data_code, is_sort1, isubcase, transient)
            self.plateStress[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strains_in_ctria3_elements(self):
        (isubcase, transient, data_code) = self.getTriHeader(True)
        data_code['table_name'] = 'OST1X'
        data = self.readTriStress(['CTRIA3'])
        if isubcase in self.plateStrain:
            self.plateStrain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            assert 'nonlinear_factor' in data_code
            self.plateStrain[isubcase] = PlateStrainObject(data_code, is_sort1, isubcase, transient)
            self.plateStrain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def getTriHeader(self, isStrain):
        """
        * analysis_code = 1 (Statics)
        * device_code   = 1 (Print)
        * table_code    = 5 (Stress)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (Real)
        * s_code        = 0 (Stress)
        * num_wide      = 8 (???)
        """
        (subcaseName, isubcase, transient, dt, analysis_code,
            is_sort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        isFiberDistance = None
        isMaxShear = None  # Von Mises/Max Shear
        if 'DISTANCE' in headers:
            isFiberDistance = True
        elif 'CURVATURE' in headers:
            isFiberDistance = False
        else:
            raise RuntimeError(headers)

        if 'VON MISES' in headers:
            isMaxShear = False
        elif 'MAX SHEAR' in headers:
            isMaxShear = True
        else:
            raise RuntimeError(headers)

        (stress_bits, s_code) = self.make_stress_bits(isFiberDistance, isMaxShear, isStrain=isStrain)
        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 5, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
                    'stress_bits': stress_bits, 'format_code': 1,
                    'element_name': 'CTRIA3', 'element_type': 74,
                    'nonlinear_factor': dt,
                    'dataNames':['lsdvmn']}
        return (isubcase, transient, data_code)

    def readTriStress(self, eType):
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
            dataTypes = [int, float, float, float, float,
                         float, float, float, float]
            sline = self.parseLine(line, dataTypes)  # line 1
            #print 'sline',sline
            sline = eType + sline
            data.append(sline)
            line = self.infile.readline()[1:].strip().split()
            #print '***',line
            dataTypes = [float, float, float, float,
                         float, float, float, float]
            sline += self.parseLine(line, dataTypes)  # line 2
            data.append(sline)
            self.i += 2
        return data

    #==========================================================================
    # CQUAD4
    def _stresses_in_cquad4_elements(self):
        (isubcase, transient, data_code) = self.getQuadHeader(2, False, 33)
        data_code['table_name'] = 'OES1X'
        data = self.readTriStress(['CQUAD4'])
        if isubcase in self.plateStress:
            self.plateStress[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            assert 'nonlinear_factor' in data_code
            self.plateStress[isubcase] = PlateStressObject(data_code, is_sort1, isubcase, transient)
            self.plateStress[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _strains_in_cquad4_elements(self):
        (isubcase, transient, data_code) = self.getQuadHeader(2, True, 33)
        data_code['table_name'] = 'OSTR1X'
        data = self.readTriStress(['CQUAD4'])
        if isubcase in self.plateStrain:
            self.plateStrain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.plateStrain[isubcase] = PlateStrainObject(data_code, is_sort1, isubcase, transient)
            self.plateStrain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def _stresses_in_cquad4_bilinear_elements(self):
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
        (isubcase, transient, data_code) = self.getQuadHeader(3, False, 144)
        #print(self.getQuadHeader(3, False, 144))
        #print("data_code =", data_code)

        data_code['table_name'] = 'OES1X'
        data = self.readQuadBilinear()
        if isubcase in self.plateStress:
            self.plateStress[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            assert 'nonlinear_factor' in data_code
            self.plateStress[isubcase] = PlateStressObject(data_code, is_sort1, isubcase, transient)
            self.plateStress[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def getQuadHeader(self, nHeaderLines, isStrain, elementNumber):
        (subcaseName, isubcase, transient, dt, analysis_code, is_sort1) = self.readSubcaseNameID()
        headers = self.skip(nHeaderLines)
        #print "headers = %s" %(headers)

        isFiberDistance = None
        isMaxShear = None  # Von Mises/Max Shear
        if 'DISTANCE' in headers:
            isFiberDistance = True
        elif 'CURVATURE' in headers:
            isFiberDistance = False
        else:
            raise RuntimeError(headers)

        if 'MAX SHEAR' in headers:
            isMaxShear = True
        elif 'VON MISES' in headers:
            isMaxShear = False
        else:
            raise RuntimeError(headers)
        (stress_bits, s_code) = self.make_stress_bits(isFiberDistance, isMaxShear, isStrain)
        data_code = {'log': self.log, 'analysis_code': analysis_code,
                    'device_code': 1, 'table_code': 5, 'sort_code': 0,
                    'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
                    'stress_bits': stress_bits, 'format_code': 1,
                    'element_name': 'CQUAD4', 'element_type': elementNumber,
                    'nonlinear_factor': dt,
                    'dataNames':['lsdvmn']}
        return (isubcase, transient, data_code)

    def readQuadBilinear(self):
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
                dataTypes = [float, float, float, float,
                             float, float, float, float]
                sline += self.parseLine(line, dataTypes)  # line 2
                data.append(sline)
                line = self.infile.readline()  # blank line
                self.i += 3

            for i in xrange(4):
                line = self.infile.readline()[1:].strip().split()
                dataTypes = [int, float, float, float, float,
                                  float, float, float, float]
                sline = self.parseLine(line, dataTypes)  # line 1
                #data.append(sline)
                line = self.infile.readline()[1:].strip().split()
                dataTypes = [float, float, float, float,
                             float, float, float, float]
                sline += self.parseLine(line, dataTypes)  # line 2
                data.append(sline)
                line = self.infile.readline()  # blank line
                self.i += 3

        return data

    #==========================================================================
    def _stresses_in_chexa_elements(self):
        return self.readSolidStress('CHEXA', 8)

    def _stresses_in_cpenta_elements(self):
        return self.readSolidStress('CPENTA', 6)

    def _stresses_in_ctetra_elements(self):
        return self.readSolidStress('CTETRA', 4)

    def _strains_in_chexa_elements(self):
        return self.readSolidStrain('CHEXA', 8)

    def _strains_in_cpenta_elements(self):
        return self.readSolidStrain('CPENTA', 6)

    def _strains_in_ctetra_elements(self):
        return self.readSolidStrain('CTETRA', 4)

    def readSolidStress(self, eType, n):
        (isubcase, transient, data_code) = self.getSolidHeader(eType, n, False)
        data_code['table_name'] = 'OES1X'

        data = self.read3DStress(eType, n)
        if isubcase in self.solidStress:
            self.solidStress[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.solidStress[isubcase] = SolidStressObject(data_code, is_sort1, isubcase, transient)
            self.solidStress[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def readSolidStrain(self, eType, n):
        (isubcase, transient, data_code) = self.getSolidHeader(eType, n, True)
        data_code['table_name'] = 'OSTR1X'

        data = self.read3DStress(eType, n)
        if isubcase in self.solidStrain:
            self.solidStrain[isubcase].add_f06_data(data, transient)
        else:
            is_sort1 = True
            self.solidStrain[isubcase] = SolidStrainObject(data_code, is_sort1, isubcase, transient)
            self.solidStrain[isubcase].add_f06_data(data, transient)
        self.iSubcases.append(isubcase)

    def getSolidHeader(self, eType, n, isStrain):
        """
        * analysis_code = 1 (Statics)
        * device_code   = 1 (Print)
        * table_code    = 5 (Stress/Strain)
        * sort_code     = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
        * format_code   = 1 (Real)
        * s_code        = 0 (Stress/Strain)
        * num_wide      = 8 (???)
        """
        (subcaseName, isubcase, transient, dt, analysis_code,
            is_sort1) = self.readSubcaseNameID()
        headers = self.skip(2)
        #print "headers = %s" %(headers)

        isMaxShear = True
        if 'VON MISES' in headers:
            isMaxShear = False

        (stress_bits, s_code) = self.make_stress_bits(
            isMaxShear=False, isStrain=isStrain, isRodOrSolid=True)
        data_code = {'log': self.log, 'analysis_code': 1, 'device_code': 1,
                    'table_code': 5, 'sort_code': 0, 'sort_bits': [0, 0, 0],
                    'num_wide': 8, 'element_name': eType, 'format_code': 1,
                    's_code': s_code, 'stress_bits': stress_bits,
                    'nonlinear_factor': dt,
                    'dataNames':['lsdvmn']}
        return (isubcase, transient, data_code)

    def read3DStress(self, eType, n):
        data = []
        while 1:
            line = self.infile.readline().rstrip('\n\r')  # [1:]
                    #              CENTER         X          #          XY             #        A         #
            sline = [line[1:17], line[17:24], line[24:28], line[28:43], line[43:47], line[47:63], line[63:66], line[66:80], line[80:83], line[83:88], line[88:93], line[93:98], line[99:113], line[113:130]]
            sline = [s.strip() for s in sline]
            if 'PAGE' in line:
                break
            elif '' is not sline[0]:
                sline = [eType] + sline
            data.append(sline)

        return data

    def make_stress_bits(self, isFiberDistance=False, isMaxShear=True, isStrain=True, isRodOrSolid=False):
        """
        ..todo:: add explanation...
        """
        #print "isMaxShear=%s isFiberDistance=%s" %(isMaxShear, isFiberDistance)

       #code = (isVonMises, isFiberCurvatur, isStress, isNotRod)
        code = (isMaxShear, isFiberDistance, isStrain, isRodOrSolid)
        mapper = {
            # element coordinate system (no material support)
            (True,  False, False,  True): ([0, 0, 0, 0, 0], 0),  # 0,  rod/csolid
            (False, False, False,  True): ([0, 0, 0, 0, 1], 1),  # 1,  rod/csolid
            (False, False,  True,  True): ([0, 0, 0, 0, 1], 1),   # ???

            (True,  False,  True, False): ([0, 1, 0, 1, 0], 10),  # 10
            (False, False,  True, False): ([0, 1, 0, 1, 1], 11),  # 11

            (True,  True,  False, False): ([0, 1, 1, 1, 0], 14),  # 14 - ???
            (True,  True,   True, False): ([0, 1, 1, 1, 0], 14),  # 14
            (False, True,   True, False): ([0, 1, 1, 1, 1], 15),  # 15

            (True, False, False, False): ([0, 0, 0, 0, 0], 0),  # 0,  composite
            (False, True, False, False): ([0, 0, 0, 0, 1], 0),  # cquad4 bilinear ??? why do i need this...
        }
        (stress_bits, s_code) = mapper[code]

        #if isMaxShear==False:
        #    stress_bits[4] = 1 # Von Mises
        #if isStrain:
        #    #stress_bits[1] = stress_bits[3] = 1 # Strain
        #    stress_bits[1] = stress_bits[3] = 1 # Strain
        #if isFiberDistance:
        #    stress_bits[2] = 1 # FiberDistance
        #print stress_bits
        #s_code = 0
        #for i,bit in enumerate(stress_bits):
        #    s_code += bit*2**i
        return (stress_bits, s_code)
