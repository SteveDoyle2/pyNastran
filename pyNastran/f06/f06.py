#pylint: disable=C0301,C0111
import os
from six.moves import zip

from pyNastran.utils import print_bad_path
from pyNastran.utils.log import get_logger

#strainEnergyDensity,TemperatureGradientObject


from pyNastran.f06.tables.oes import OES
from pyNastran.f06.tables.oug import OUG
from pyNastran.f06.tables.oqg import OQG
from pyNastran.f06.tables.oef import OEF
from pyNastran.f06.tables.lama import LAMA
from pyNastran.f06.tables.max_min import MAX_MIN
from pyNastran.f06.f06Writer import F06Writer
from pyNastran.op2.tables.ogf_gridPointForces.ogf_Objects import RealGridPointForces

from pyNastran.utils import is_binary
from pyNastran.f06.errors import FatalError

#class F06Deprecated(object):
    #def __init__(self, f06_filename):
        #self.f06FileName = f06_filename
        #self.f06_filename = None
    #def read_f06(self, f06_filename):
        #pass
    #def readF06(self):
        #"""... seealso::: read_f06"""
        #self.read_f06(self.f06_filename)


class F06(OES, OEF, OUG, OQG, LAMA, MAX_MIN, F06Writer): #, F06Deprecated):
    def stop_after_reading_grid_point_weight(self, stop=True):
        self._stop_after_reading_mass = stop

    def __init__(self, debug=False, log=None):
        """
        Initializes the F06 object

        :makeGeom:    reads the BDF tables (default=False)
        :debug:       prints data about how the F06 was parsed (default=False)
        :log:         a logging object to write debug messages to

        .. seealso:: import logging
        """
        OES.__init__(self)
        OEF.__init__(self)
        OQG.__init__(self)
        OUG.__init__(self)
        LAMA.__init__(self)
        MAX_MIN.__init__(self)
        #F06Deprecated.__init__(self, f06_filename)
        #F06Writer.__init__(self)

        self._subtitle = None
        self.card_count = {}
        self._stop_after_reading_mass = False
        self.stored_lines = []
        self.i = 0

        self.__init_data__(debug, log)

        self._line_marker_map = {
            'R E A L   E I G E N V E C T O R   N O' : self._real_eigenvectors,
            'C O M P L E X   E I G E N V E C T O R   NO' : self._complex_eigenvectors,
            'News file -' : self._executive_control_echo,
        }
        self._marker_map = {
            #====================================================================
            # debug info
            'N A S T R A N    F I L E    A N D    S Y S T E M    P A R A M E T E R    E C H O' : self._nastran_file_and_system_parameter_echo,
            'N A S T R A N    E X E C U T I V E    C O N T R O L    E C H O' : self._executive_control_echo,
            'C A S E    C O N T R O L    E C H O' : self._case_control_echo,
            'G R I D   P O I N T   S I N G U L A R I T Y   T A B L E' : self._grid_point_singularity_table,

            # dummy
            'E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y' : self._executive_control_echo,
            'M O D E L   S U M M A R Y' : self._executive_control_echo,
            'M O D E L   S U M M A R Y          BULK = 0' : self._executive_control_echo,
            'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._grid_point_singularity_table,
            'G R I D   P O I N T   F O R C E   B A L A N C E' : self._grid_point_force_balance,
            'N A S T R A N   S O U R C E   P R O G R A M   C O M P I L A T I O N             SUBDMAP  =  SESTATIC' : self._executive_control_echo,
            #====================================================================
            # useful info
            #'E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y'
            'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R': self._grid_point_weight_generator,

            # dummy
            'MAXIMUM  SPCFORCES' : self._get_max_spc_forces,
            #'OLOAD    RESULTANT' : self.getMaxMpcForces,
            'MAXIMUM  MPCFORCES' : self._get_max_mpc_forces,
            'SPCFORCE RESULTANT' : self._get_max_mpc_forces,
            'MPCFORCE RESULTANT' : self._get_max_mpc_forces,
            'MAXIMUM  DISPLACEMENTS' : self._get_max_displacements,
            'MAXIMUM  APPLIED LOADS' : self._get_max_applied_loads,


            #====================================================================
            # F06 specific tables
            #'N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S' : self._nondimensional_stability_and_control_deriv_coeffs,
            #'N O N - D I M E N S I O N A L    H I N G E    M O M E N T    D E R I V A T I V E   C O E F F I C I E N T S':  self._nondimensional_hinge_moment_derivative_coeffs,
            #'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S': self._aerostatic_data_recovery_output_tables,
            #'S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S': self._structural_monitor_point_integrated_loads,

            'N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S' : self._executive_control_echo,
            'N O N - D I M E N S I O N A L    H I N G E    M O M E N T    D E R I V A T I V E   C O E F F I C I E N T S':  self._executive_control_echo,
            'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S': self._executive_control_echo,
            'S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S': self._executive_control_echo,
            'FLUTTER  SUMMARY' : self._executive_control_echo,
            #------------------------
            #'R O T O R   D Y N A M I C S   S U M M A R Y'
            #'R O T O R   D Y N A M I C S   M A S S   S U M M A R Y'
            #'E I G E N V A L U E  A N A L Y S I S   S U M M A R Y   (COMPLEX LANCZOS METHOD)'

            'R O T O R   D Y N A M I C S   S U M M A R Y': self._executive_control_echo,
            'R O T O R   D Y N A M I C S   M A S S   S U M M A R Y': self._executive_control_echo,
            'E I G E N V A L U E  A N A L Y S I S   S U M M A R Y   (COMPLEX LANCZOS METHOD)': self._executive_control_echo,

            #------------------------
            #====================================================================

            # OUG tables
            'R E A L   E I G E N V A L U E S': self._real_eigenvalues,
            'C O M P L E X   E I G E N V A L U E   S U M M A R Y':self._complex_eigenvalue_summary,

            'D I S P L A C E M E N T   V E C T O R' : self._displacement_vector,
            'C O M P L E X   D I S P L A C E M E N T   V E C T O R' : self._complex_displacement_vector,
            'F O R C E S   O F   S I N G L E - P O I N T   C O N S T R A I N T' : self._forces_of_single_point_constraints,
            'F O R C E S   O F   M U L T I P O I N T   C O N S T R A I N T' : self._forces_of_multi_point_constraints,

            'T E M P E R A T U R E   V E C T O R' : self._temperature_vector,
            'L O A D   V E C T O R' : self._load_vector,
            #'F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S': self._temperature_gradients_and_fluxes,

            #====================================================================
            # OES O-D
            'S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )' : self._stresses_in_cbar_elements,
            'S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )'  : self._strains_in_cbar_elements,

            'S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )' : self._stresses_in_crod_elements,
            'S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )'   : self._strains_in_crod_elements,

            'S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )' : self._stresses_in_crod_elements,
            'S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )'   : self._strains_in_crod_elements,
            #====================================================================
            # OES 1-D
            'S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 1 )' : self._stresses_in_celas1_elements,
            'S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )' : self._stresses_in_celas2_elements,
            'S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 3 )' : self._stresses_in_celas3_elements,
            'S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 4 )' : self._stresses_in_celas4_elements,

            'S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 1 )' : self._strains_in_celas1_elements,
            'S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 2 )' : self._strains_in_celas2_elements,
            'S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 3 )' : self._strains_in_celas3_elements,
            'S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 4 )' : self._strains_in_celas4_elements,
            #====================================================================
            # OES 2-D (no support for CQUAD8, CQUAD, CQUADR, CTRIAR, CTRAI6, CTRIAX, CTRIAX6)
            'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._stresses_in_ctria3_elements,
            'S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._strains_in_ctria3_elements,

            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._stresses_in_cquad4_elements,
            'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._strains_in_cquad4_elements,

            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN' : self._stresses_in_cquad4_bilinear_elements,
            'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN' : self._strains_in_cquad4_bilinear_elements,

            #==
            # composite partial ??? (e.g. bilinear quad)
            'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )': self._stresses_in_composite_ctria3_elements,
            'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )': self._stresses_in_composite_cquad4_elements,

            'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )' : self._strains_in_composite_ctria3_elements,
            'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )' : self._strains_in_composite_cquad4_elements,

            #===
            # ..todo:: not implemented
            'S T R E S S E S   I N   T R I A X 6   E L E M E N T S' : self._executive_control_echo,
            #====================================================================
            # OES 3-D
            'S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )' : self._stresses_in_ctetra_elements,
            'S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._stresses_in_chexa_elements,
            'S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._stresses_in_cpenta_elements,

            'S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )' : self._strains_in_ctetra_elements,
            'S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._strains_in_chexa_elements,
            'S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._strains_in_cpenta_elements,
            #====================================================================
            # more not implemented...

            # STRESS
            'S T R E S S E S   I N   H Y P E R E L A S T I C   H E X A H E D R O N   E L E M E N T S  ( HEXA8FD )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S    ( Q U A D 4 )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   Q U A D R I L A T E R A L   E L E M E N T S  ( QUAD4FD )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S  IN  H Y P E R E L A S T I C   A X I S Y M M.  Q U A D R I L A T E R A L  ELEMENTS (QUADXFD)' : self._executive_control_echo,

            # FORCE
            'F O R C E S   I N   R O D   E L E M E N T S     ( C R O D )': self._forces_in_crod_elements,
            'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN':  self._forces_in_cquad4s_bilinear,

            'F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )' : self._executive_control_echo,
            'F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )':  self._executive_control_echo,
            'F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 1 )': self._executive_control_echo,
            'F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )': self._executive_control_echo,
            'F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 3 )': self._executive_control_echo,
            'F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 4 )': self._executive_control_echo,

            # energy
            'E L E M E N T   S T R A I N   E N E R G I E S' : self._element_strain_energies,
            'E L E M E N T   S T R A I N   E N E R G I E S   ( A V E R A G E )' : self._executive_control_echo,

            # complex stress
            'C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )' : self._executive_control_echo,
            #'* * * END OF JOB * * *': self.end(),
        }
        self.markers = self._marker_map.keys()

    def __init_data__(self, debug=False, log=None):
        OES.__init__(self)
        OQG.__init__(self)
        OUG.__init__(self)
        F06Writer.__init__(self)

        #: the TITLE in the Case Control Deck
        self.Title = ''
        self._start_log(log, debug)

    def _start_log(self, log=None, debug=False):
        """
        Sets up a dummy logger if one is not provided

        :self:  the object pointer
        :log:   a python logging object
        :debug: adds debug messages (True/False)
        """
        self.log = get_logger(log, 'debug' if debug else 'info')

    def _get_grid_point_singularities(self):  # .. todo:: not done
        """
        ::

                      G R I D   P O I N T   S I N G U L A R I T Y   T A B L E
          POINT    TYPE   FAILED      STIFFNESS       OLD USET           NEW USET
           ID            DIRECTION      RATIO     EXCLUSIVE  UNION   EXCLUSIVE  UNION
            1        G      4         0.00E+00          B        F         SB       S    *
            1        G      5         0.00E+00          B        F         SB       S    *
        """
        pass

    def _grid_point_weight_generator(self):
        line = ''
        lines = []
        while 'PAGE' not in line:
            line = self.infile.readline()[1:].strip()
            lines.append(line)
            self.i += 1
            self.fatal_check(line)
        #print('\n'.join(lines))
        self.grid_point_weight.read_grid_point_weight(lines)

    def _case_control_echo(self):
        line = ''
        lines = []
        while 'PAGE' not in line:
            line = self.infile.readline()[1:].strip()
            lines.append(line)
            self.i += 1
            self.fatal_check(line)
        #self.grid_point_weight.read_grid_point_weight(lines)

    def _executive_control_echo(self):
        line = ''
        lines = []
        while 'PAGE' not in line:
            line = self.infile.readline()[1:].strip()
            lines.append(line)
            self.i += 1
            self.fatal_check(line)
        #self.grid_point_weight.read_grid_point_weight(lines)

    def fatal_check(self, line):
        if 'FATAL' in line:
            raise FatalError(line)

    def _nastran_file_and_system_parameter_echo(self):
        line = ''
        lines = []
        while 'PAGE' not in line:
            line = self.infile.readline()[1:].strip()
            lines.append(line)
            self.i += 1
            self.fatal_check(line)
        #self.grid_point_weight.read_grid_point_weight(lines)

    def _grid_point_singularity_table(self):
        line = ''
        lines = []
        while 'PAGE' not in line:
            line = self.infile.readline()[1:].strip()
            lines.append(line)
            self.i += 1
            self.fatal_check(line)
        #self.grid_point_weight.read_grid_point_weight(lines)

    def _set_f06_date(self, month, day, year):
        months = ['JANUARY', 'FEBRUARY', 'MARCH', 'APRIL', 'MAY', 'JUNE',
                  'JULY', 'AUGUST', 'SEPTEMBER', 'OCTOBER', 'NOVEMBER', 'DECEMBER']
        assert month in months, 'month=%r' % month
        month = months.index(month) + 1
        day = int(day)
        year = int(year)
        self.date = (month, day, year)

    def _read_f06_subcase_header(self, n=-2):
        """
        -4 -> 1                                                     JANUARY   5, 2014  MSC.NASTRAN 11/25/11   PAGE    14
        -3 -> DEFAULT
        -2 -> xxx             subcase 1
        """
        #print('-----------------')
        ### changed...
        #for iline in [-4, -3, -2, -1, 0]:
            #if iline == -1:
                #print('***line[%s]=%r' % (iline, self.stored_lines[iline].replace("    ", " ")))
            #else:
                #print('line[%s]=%r' % (iline, self.stored_lines[iline].replace("    ", " ")))
        #assert 'PAGE' in self.stored_lines[-4], self.stored_lines[-4]
        subtitle = self.stored_lines[-3].strip()

        msg = ''
        lines2 = []
        for i, line in enumerate(self.stored_lines[-4:]):
            line2 = line.rstrip()
            if line2:
                #msg += '%i -> %s\n' % (-4 + i, line.rstrip())
                #print("%r" % line2.replace("  ", " "))
                lines2.append(line2)

        if self.Title is None or self.Title == '' and len(self.stored_lines) > 4:
            title_line = self.stored_lines[-4]
            self.Title = title_line[1:75].strip()
            date = title_line[75:93].strip()
            if date:
                try:
                    month, day, year = date.split()
                except:
                    raise RuntimeError('Couldnt parse date; line=%r' % title_line.strip())
                self._set_f06_date(month, day[:-1], year)  # -1 chops the comma
            #assert 'PAGE' not in title_line, '%r' % date
            assert 'D I S P L A C' not in self.Title, msg
        #self.Title = subcaseName  # 'no title'

        subcase_name = ''
        #print("subcaseLine = %r" % subcaseName)

        # hack...
        if 'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN' in self.stored_lines[0]:
            label_isubcase = self.stored_lines[-3]
        else:
            label_isubcase = self.stored_lines[-2]
        label, isubcase = _parse_label_isubcase(label_isubcase)

        #subtitle = 'SUBCASE %s' % isubcase
        #label = 'SUBCASE %s' % isubcase

        #self.iSubcaseNameMap[self.isubcase] = [self.subtitle, self.label]

#title      date_stamp  page_stamp
#subtitle
#label      ???

        #if 'F O R C E S' in self.stored_lines[0]:
            #print("title    = %r" % self.Title)
            #print("subtitle = %r" % subtitle)
            #print("label    = %r" % label)
            #print("isubcase = %r" % isubcase)
        #print('------------')

        #assert self.Title == 'MSC.NASTRAN JOB CREATED ON 12-MAR-13 AT 12:52:23', self.Title
        self._subtitle = subtitle
        self.iSubcaseNameMap[isubcase] = [subtitle, label]
        transient = self.stored_lines[-1].strip()
        is_sort1 = False
        if transient:
            try:
                trans_word, trans_value = transient.split('=')
            except ValueError:
                msg = 'transient = %r' % transient
                msg += 'stored lines = [%r]' % '\n'.join(self.stored_lines)
                #print(msg)
                raise ValueError(msg)
            trans_word = trans_word.strip()
            trans_value = float(trans_value)
            transient = [trans_word, trans_value]

            if trans_word == 'LOAD STEP':  # nonlinear statics
                analysis_code = 10
            elif trans_word == 'TIME STEP':  # TODO check name
                analysis_code = 6
            elif trans_word == 'EIGENVALUE':  # normal modes
                analysis_code = 2
            elif trans_word == 'FREQ':  # TODO check name
                analysis_code = 5
            elif trans_word == 'FREQUENCY':
                analysis_code = 5
            elif trans_word == 'POINT-ID':
                is_sort1 = True
                analysis_code = None
            elif trans_word == 'ELEMENT-ID':
                is_sort1 = True
                #is_sort1 = False
                #is_sort2 = True
                analysis_code = None
            else:
                raise NotImplementedError('transientWord=%r is not supported...' % trans_word)
        else:
            transient = None
            analysis_code = 1

        dt = None
        if transient is not None:
            dt = transient[1]
        return (subcase_name, isubcase, transient, dt, analysis_code, is_sort1)

    def _element_strain_energies(self):
        """
        ::

          EIGENVALUE = -3.741384E-04
          CYCLES =  3.078479E-03
                                             E L E M E N T   S T R A I N   E N E R G I E S

                  ELEMENT-TYPE = QUAD4               * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM     =  -1.188367E-05
                     MODE               1            * TOTAL ENERGY OF ALL ELEMENTS IN SET      -1 =  -1.188367E-05

                                      ELEMENT-ID          STRAIN-ENERGY           PERCENT OF TOTAL    STRAIN-ENERGY-DENSITY
                                               1         -5.410134E-08                -0.0929             -4.328107E-05
                                               2         -3.301516E-09                -0.0057             -2.641213E-06
        """
        isubcase = 1 # TODO not correct
        cycles = self.stored_lines[-1][1:].strip()
        try:
            cycles = float(cycles.split('=')[1])
        except (IndexError, ValueError):
            return

        eigenvalue = self.stored_lines[-2][1:].strip()
        eigenvalue = float(eigenvalue.split('=')[1])
        #print("eigenvalue=%s cycle=%s" % (eigenvalue, cycles))

        eTypeLine = self.skip(2)[1:]
        eType = eTypeLine[30:40]
        totalEnergy1 = eTypeLine[99:114]

        mode_line = self.skip(1)[1:]
        iMode = mode_line[24:40]
        totalEnergy2 = mode_line[99:114]
        #print("eType=%r totalEnergy1=%r" % (eType, totalEnergy1))
        #print("iMode=%r totalEnergy2=%r" % (iMode, totalEnergy2))
        headers = self.skip(2)

        data = []
        while 1:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            self.i += 1
            if 'PAGE' in line:
                break
            sline = line.strip().split()
            if sline == []:
                break
            eid = int(sline[0])
            strain_energy = float(sline[1])
            percent_total = float(sline[2])
            strain_energy_density = float(sline[3])
            out = (eid, strain_energy, percent_total, strain_energy_density)
            data.append(out)

        if sline == []:
            line = self.infile.readline()[1:].rstrip('\r\n ')
            self.i += 1

        return
        if isubcase in self.iSubcases:
            self.strainEnergyDensity[isubcase].readF06Data(data, transient)
        else:
            sed = strain_energy_density(data, transient)
            sed.readF06Data(data, transient)
            self.strainEnergyDensity[isubcase] = sed

    def _grid_point_force_balance(self):
        try:
            (subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
        except ValueError:
            line = ''
            while 'PAGE' not in line:
                line = self.infile.readline()[1:].rstrip()
                #lines.append(line)
                self.i += 1
            return

        headers = self.skip(2)
        line = ''
        lines = []
        #POINT-ID ELEMENT-ID SOURCE T1 T2 T3 R1 R2 R3
        while 'PAGE' not in line:
            line = self.infile.readline()[1:].rstrip()
            lines.append(line)
            self.i += 1
            self.fatal_check(line)

        data = []
        for line in lines:
            if 'PAGE' in line:
                break
            point_id, element_id, source = line[1:11], line[11:24], line[24:36]
            t1, t2, t3, r1, r2, r3 = line[36:55], line[55:72], line[72:87], line[87:102], line[102:117], line[117:138]
            try:
                point_id = int(point_id)
                element_id = element_id.strip()
                if element_id:
                    element_id = int(element_id)
                t1 = float(t1)
                t2 = float(t2)
                t3 = float(t3)
                r1 = float(r1)
                r2 = float(r2)
                r3 = float(r3)
            except:
                msg = "point_id=%r, element_id=%r, source=%r\n" % (point_id, element_id, source)
                msg += "t1=%r, t2=%r, t3=%r, r1=%r, r2=%r, r3=%r" % (t1, t2, t3, r1, r2, r3)
                raise SyntaxError(msg)

            data.append([point_id, element_id, source, t1, t2, t3, r1, r2, r3])

            data_code = {'analysis_code': analysis_code,
                        'device_code': 1, 'sort_code': 0,
                        'sort_bits': [0, 0, 0], 'num_wide': 9,

                        'table_code': 1,
                        'table_name': 'OES1X',

                        #'s_code': s_code,
                        #'stress_bits': stress_bits,
                        #'element_name': element_name, 'element_type': element_type,

                        'format_code': 1,
                        'nonlinear_factor': dt,
                        'dataNames':['lsdvmn'],
                        'lsdvmn': 1,
                        }
        is_sort1 = True
        if isubcase not in self.gridPointForces:
            self.gridPointForces[isubcase] = RealGridPointForces(data_code, is_sort1, isubcase, dt)
        self.gridPointForces[isubcase].add_f06_data(dt, data)
        self.iSubcases.append(isubcase)

    def _read_f06_table(self, Format, debug=False):
        """
        Reads displacement, spc/mpc forces

        :self:   the object pointer
        :Format: .. seealso:: parseLine
        """
        sline = True
        data = []
        while sline:
            line = self.infile.readline()[1:].strip()
            sline = line.split()
            if debug:
                print(sline)
            self.i += 1
            if 'PAGE' in sline:
                #self.stored_lines = [line]  ## changed...
                return data
            sline = self.parseLine(sline, Format)
            if sline is None or len(sline) == 0:
                return data
            data.append(sline)
        return data

    def _read_table_dummy(self):
        sline = True
        data = []
        while sline:
            sline = self.infile.readline()[1:].strip().split()
            self.i += 1
            if 'PAGE' in sline:
                return data
            if sline is None:
                return data
            data.append(sline)
        return data

    def parseLine(self, sline, Format):
        """
        :self:   the object pointer
        :sline:  list of strings (split line)
        :Format: list of types [int,str,float,float,float] that maps to sline
        """
        out = []
        for entry, iFormat in zip(sline, Format):
            try:
                entry2 = iFormat(entry)
            except:
                #print("sline=%r\n entry=%r format=%s" % (sline, entry, iFormat))
                #raise
                return None
            out.append(entry2)
        return out

    def _parse_line_blanks(self, sline, Format):
        """allows blanks"""
        out = []

        for entry, iFormat in zip(sline, Format):
            if entry.strip():
                try:
                    entry2 = iFormat(entry)
                except:
                    print("sline=%r\n entry=%r format=%s" % (sline, entry, Format))
                    raise
            else:
                entry2 = None
                #print("sline=%r\n entry=%r format=%s" % (sline, entry, iFormat))
            out.append(entry2)
        return out

    def read_f06(self, f06_filename=None):
        """
        Reads the F06 file

        :self: the object pointer

        :f06FileName: the file to be parsed (None -> GUI)
        """
        if f06_filename is None:
            from pyNastran.utils.gui_io import load_file_dialog
            wildcard_wx = "Nastran F06 (*.f06)|*.f06|" \
                "All files (*.*)|*.*"
            wildcard_qt = "Nastran F06 (*.f06);;All files (*)"
            title = 'Please select a F06 to load'
            f06_filename = load_file_dialog(title, wildcard_wx, wildcard_qt)
            assert f06_filename is not None, f06_filename
        else:
            if not os.path.exists(f06_filename):
                msg = 'cant find f06_filename=%r\n%s' \
                    % (f06_filename, print_bad_path(f06_filename))
                raise RuntimeError(msg)

        if is_binary(f06_filename):
            raise IOError('f06_filename=%r is not a binary F06.' % f06_filename)
        if os.path.getsize(f06_filename) == 0:
            raise IOError('f06_filename=%r is empty.' % f06_filename)

        self.log.debug('f06_filename = %r' % f06_filename)
        self.f06_filename = f06_filename


        blank = 0
        self.infile = open(self.f06_filename, 'r')
        while 1:
            #if self.i % 1000 == 0:
                #print("i=%i" % (self.i))
            line = self.infile.readline()
            marker = line[1:].strip()

            if 'FATAL' in marker and 'IF THE FLAG IS FATAL' not in marker:
                msg = '\n' + marker
                fatal_count = 0
                while 1:
                    line = self.infile.readline().rstrip()
                    #print("blank = %s" % blank)
                    fatal_count += 1
                    if fatal_count == 20 or '* * * END OF JOB * * *' in line:
                        break
                    #else:
                        #blank = 0
                    msg += line + '\n'
                raise FatalError(msg.rstrip())

            if(marker != '' and 'SUBCASE' not in marker and 'PAGE' not in marker and 'FORTRAN' not in marker
               and 'USER INFORMATION MESSAGE' not in marker and 'TOTAL DATA WRITTEN FOR DATA BLOCK' not in marker
               and marker not in self.markers and marker != self._subtitle):
                #print("marker = %r" % marker)
                pass
                #print('Title  = %r' % self.subtitle)

            if marker in self.markers:
                blank = 0
                #print("\n1*marker = %r" % marker)
                self._marker_map[marker]()
                if(self._stop_after_reading_mass and
                   marker in 'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R'):
                    break
                self.stored_lines = []
            elif 'R E A L   E I G E N V E C T O R   N O' in marker:
                blank = 0
                #print("\n2*marker = %r" % marker)
                self._line_marker_map['R E A L   E I G E N V E C T O R   N O'](marker)
                #print('done with real eigenvector')
                self.stored_lines = []
            elif 'C O M P L E X   E I G E N V E C T O R   NO' in marker:
                blank = 0
                #print("\n2*marker = %r" % marker)
                self._line_marker_map['C O M P L E X   E I G E N V E C T O R   NO'](marker)
                self.stored_lines = []

            elif 'News file -' in marker:
                blank = 0
                self._line_marker_map['News file -']()
                self.stored_lines = []
            elif marker == '':
                blank += 1
                if blank == 20:
                    break
            elif self._is_marker(marker):  # marker with space in it (e.g. Model Summary)
                print("***marker = %r" % marker)

            else:
                blank = 0

            self.stored_lines.append(line)
            self.i += 1
        #print("i=%i" % self.i)
        self.infile.close()
        self._process_f06()
        if hasattr(self, '_ieigenvalue'):
            del self._ieigenvalue

    def _process_f06(self):
        #data = [self.disp,self.SpcForces,self.stress,self.isoStress,self.barStress,self.solidStress,self.temperature]
        data_pack = [self.solidStress]
        for data_set in data_pack:
            for key, data in data_set.iteritems():
                data.processF06Data()

    def _is_marker(self, marker):
        """returns True if the word follows the 'N A S T R A N   P A T T E R N'"""
        marker = marker.strip().split('$')[0].strip()

        if len(marker) < 2 or marker == '* * * * * * * * * * * * * * * * * * * *':
            return False
        for i, char in enumerate(marker):
            #print("i=%s i%%2=%s char=%s" % (i, i%2, char))
            if i % 2 == 1 and ' ' is not char:
                return False
            elif i % 2 == 0 and ' ' == char:
                return False
        return True

    def skip(self, iskip):
        for i in xrange(iskip - 1):
            self.infile.readline()
        self.i += iskip
        return self.infile.readline()

    #def get_op2_stats(self):
        #"""
        #Gets info about the contents of the different attributes of the
        #OP2 class.
        #"""
        #pass

def _parse_label_isubcase(label_isubcase):
    label = label_isubcase[1:65].strip()
    isubcase = label_isubcase[65:].strip()
    #print('stored2 = ', label_isubcase.strip().replace('         ', ' '))
    if isubcase:
        isubcase = int(isubcase.split()[-1])
    else:
        #raise RuntimeError('asdf')
        isubcase = 1
    return label, isubcase
    #assert isinstance(isubcase,int),'isubcase=%r' % isubcase
    #print("subcaseName=%s isubcase=%s" % (subcaseName, isubcase))

if __name__ == '__main__':
    from pyNastran.f06.test.test_f06 import main
    main()

