#pylint: disable=C0301,C0111
from __future__ import print_function
from six import iteritems
from six.moves import zip, range
import os
from codecs import open

from pyNastran import is_release
from pyNastran.utils import print_bad_path
from pyNastran.utils.log import get_logger

#strainEnergyDensity,TemperatureGradientObject
from pyNastran.op2.tables.oee_energy.oee_objects import RealStrainEnergy

from pyNastran.f06.tables.oes import OES
from pyNastran.f06.tables.oug import OUG
from pyNastran.f06.tables.oqg import OQG
from pyNastran.f06.tables.oef import OEF
from pyNastran.f06.tables.lama import LAMA
from pyNastran.f06.tables.max_min import MAX_MIN
from pyNastran.f06.f06_writer import F06Writer
from pyNastran.op2.tables.ogf_gridPointForces.ogf_Objects import RealGridPointForces

from pyNastran.utils import is_binary_file
from pyNastran.f06.errors import FatalError


class F06(OES, OEF, OUG, OQG, LAMA, MAX_MIN, F06Writer):
    def stop_after_reading_grid_point_weight(self, stop=True):
        self._stop_after_reading_mass = stop

    def build_vectorization(self):
        if self.is_vectorized:
            table_types = self.get_table_types()
            for table_type in table_types:
                result = getattr(self, table_type)
                for key, case in iteritems(result):
                    if hasattr(case, 'build_f06_vectorization'):
                        case.build_f06_vectorization()

    def __init__(self, debug=False, log=None):
        """
        Initializes the F06 object

        :param makeGeom:    reads the BDF tables (default=False)
        :param debug:       prints data about how the F06 was parsed (default=False)
        :param log:         a logging object to write debug messages to

        .. seealso:: import logging
        """
        OES.__init__(self)
        OEF.__init__(self)
        OQG.__init__(self)
        OUG.__init__(self)
        LAMA.__init__(self)
        MAX_MIN.__init__(self)
        F06Writer.__init__(self)

        self.f06_filename = None
        self._subtitle = None
        self.card_count = {}
        self._stop_after_reading_mass = False
        self.stored_lines = []
        self.i = 0

        #: the TITLE in the Case Control Deck
        self.title = ''
        self._start_log(log, debug)

        self._line_marker_map = {
            'R E A L   E I G E N V E C T O R   N O' : self._real_eigenvectors,
            'C O M P L E X   E I G E N V E C T O R   NO' : self._complex_eigenvectors,
            'News file -' : self._executive_control_echo,
        }
        self._marker_map = {
            #====================================================================
            # debug info
            'THIS PROGRAM IS CONFIDENTIAL AND A TRADE SECRET OF MSC.SOFTWARE CORPORATION.  THE RECEIPT OR' : self._executive_control_echo,
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

            'U S E T   D E F I N I T I O N   T A B L E   ( I N T E R N A L   S E Q U E N C E ,   C O L U M N   S O R T )' : self._executive_control_echo,
            'N O N - L I N E A R   I T E R A T I O N   M O D U L E   O U T P U T' : self._executive_control_echo,
            'N A S T R A N   D M A P / N D D L   L I N K A G E   E D I T O R' : self._executive_control_echo,
            'N O N - L I N E A R   I T E R A T I O N   M O D U L E   S O L U T I O N   C O N T R O L   D A T A' : self._executive_control_echo,
            #====================================================================
            # useful info
            #'E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y'
            'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R': self._grid_point_weight_generator,

            # dummy
            'MAXIMUM  SPCFORCES' : self._get_max_spc_forces,
            'OLOAD    RESULTANT' : self._get_oload_resultant,
            'MAXIMUM  MPCFORCES' : self._get_max_mpc_forces,
            'SPCFORCE RESULTANT' : self._get_max_mpc_forces,
            'MPCFORCE RESULTANT' : self._get_max_mpc_forces,
            'MAXIMUM  DISPLACEMENTS' : self._get_max_displacements,
            'MAXIMUM  APPLIED LOADS' : self._get_max_applied_loads,


            #====================================================================
            # F06 specific tables
            'S O R T E D   B U L K   D A T A   E C H O' : self._executive_control_echo,
            'I N P U T   B U L K   D A T A   E C H O' : self._executive_control_echo,
            #'N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S' : self._nondimensional_stability_and_control_deriv_coeffs,
            #'N O N - D I M E N S I O N A L    H I N G E    M O M E N T    D E R I V A T I V E   C O E F F I C I E N T S' : self._nondimensional_hinge_moment_derivative_coeffs,
            #'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S' : self._aerostatic_data_recovery_output_tables,
            #'S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S' : self._structural_monitor_point_integrated_loads,

            'N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S' : self._executive_control_echo,
            'N O N - D I M E N S I O N A L    H I N G E    M O M E N T    D E R I V A T I V E   C O E F F I C I E N T S' :  self._executive_control_echo,
            'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S' : self._executive_control_echo,
            'S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S' : self._executive_control_echo,

            #'A E R O D Y N A M I C   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S' : self._executive_control_echo,  # uncommmented => crash
            'FLUTTER  SUMMARY' : self._executive_control_echo,
            #------------------------
            #'R O T O R   D Y N A M I C S   S U M M A R Y' : self._executive_control_echo,
            #'R O T O R   D Y N A M I C S   M A S S   S U M M A R Y' : self._executive_control_echo,


            #'R O T O R   D Y N A M I C S   S U M M A R Y' : self._executive_control_echo,  # uncommmented => crash
            #'R O T O R   D Y N A M I C S   M A S S   S U M M A R Y' : self._executive_control_echo,  # uncommmented => crash
            'E I G E N V A L U E  A N A L Y S I S   S U M M A R Y   (COMPLEX LANCZOS METHOD)' : self._executive_control_echo,
            #'E I G E N V A L U E  A N A L Y S I S   S U M M A R Y   (READ MODULE)' : self._executive_control_echo,

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


            # stress - good
            'S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 1 )' : self._stress_in_celas1_elements,
            'S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )' : self._stress_in_celas2_elements,
            'S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 3 )' : self._stress_in_celas3_elements,
            'S T R E S S E S   I N   S C A L A R   S P R I N G S        ( C E L A S 4 )' : self._stress_in_celas4_elements,

            'S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )' : self._stress_in_crod_elements,
            'S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )' : self._stress_in_ctube_elements,
            'S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )' : self._stress_in_conrod_elements,
            'S T R E S S E S   I N   B A R   E L E M E N T S          ( C B A R )' : self._stress_in_cbar_elements,

            'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._stress_in_ctria3_elements,
            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._stress_in_cquad4_elements,
            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN' : self._stress_in_cquad4_bilinear_elements,

            'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )': self._stress_in_composite_ctria3_elements,
            'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )': self._stress_in_composite_cquad4_elements,

            'S T R E S S E S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )' : self._stress_in_ctetra_elements,
            'S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._stress_in_chexa_elements,
            'S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._stress_in_cpenta_elements,

            # strain - good
            'S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 1 )' : self._strain_in_celas1_elements,
            'S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 2 )' : self._strain_in_celas2_elements,
            'S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 3 )' : self._strain_in_celas3_elements,
            'S T R A I N S    I N   S C A L A R   S P R I N G S        ( C E L A S 4 )' : self._strain_in_celas4_elements,

            'S T R A I N S   I N   R O D   E L E M E N T S      ( C R O D )' : self._strain_in_crod_elements,
            'S T R A I N S   I N   R O D   E L E M E N T S      ( C O N R O D )' : self._strain_in_conrod_elements,
            'S T R A I N S    I N   B A R   E L E M E N T S          ( C B A R )' : self._strain_in_cbar_elements,

            'S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._strain_in_ctria3_elements,
            'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._strain_in_cquad4_elements,
            'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN' : self._strain_in_cquad4_bilinear_elements,

            'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )' : self._strain_in_composite_ctria3_elements,
            'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )' : self._strain_in_composite_cquad4_elements,

            'S T R A I N S   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )' : self._strain_in_ctetra_elements,
            'S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._strain_in_chexa_elements,
            'S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._strain_in_cpenta_elements,

            # energy - good
            'E L E M E N T   S T R A I N   E N E R G I E S' : self._element_strain_energies,

            #====================================================================

            # stress - not implemented
            'S T R E S S   D I S T R I B U T I O N   I N   B A R   E L E M E N T S       ( C B A R )' : self._executive_control_echo,
            'S T R E S S E S   I N   B E N D   E L E M E N T S        ( C B E N D )' : self._executive_control_echo,
            'S T R E S S E S   I N   B E A M   E L E M E N T S        ( C B E A M )' : self._executive_control_echo,

            'S T R E S S E S   ( F O R C E S )   I N   B U S H 1 D   E L E M E N T S   ( C B U S H 1 D )' : self._executive_control_echo,
            'S T R E S S E S   I N   B U S H   E L E M E N T S        ( C B U S H )' : self._executive_control_echo,

            'S T R E S S E S   I N   T R I A X 6   E L E M E N T S' : self._executive_control_echo,
            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
            'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
            'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,

            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = CENTER' : self._executive_control_echo,
            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = CUBIC' : self._executive_control_echo,
            'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = SGAGE' : self._executive_control_echo,

            'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
            'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,
            'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
            'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,

            # strain - not implemented
            'S T R A I N S   I N   B U S H   E L E M E N T S        ( C B U S H )' : self._executive_control_echo,

            'S T R A I N   D I S T R I B U T I O N   I N   B A R   E L E M E N T S       ( C B A R )' : self._executive_control_echo,
            'S T R A I N S    I N   B E A M   E L E M E N T S        ( C B E A M )' : self._executive_control_echo,
            'S T R A I N S    I N   B E N D   E L E M E N T S        ( C B E N D )' : self._executive_control_echo,

            'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
            'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
            'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = CENTER' : self._executive_control_echo,
            'S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
            'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,

            # other real - not implemented
            'S T R E S S E S   I N   U S E R   E L E M E N T S (CDUM8)' : self._executive_control_echo,
            'S T R E S S E S   I N   U S E R   E L E M E N T S (CDUM9)' : self._executive_control_echo,
            'E L E M E N T   S T R A I N   E N E R G I E S   ( A V E R A G E )' : self._executive_control_echo,
            'E L E M E N T   S T R E S S   D I S C O N T I N U I T I E S  - -     S U R F A C E      99' : self._executive_control_echo,

            # FORCE - not implemented
            'F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 1 )' : self._executive_control_echo,
            'F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 2 )' : self._executive_control_echo,
            'F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 3 )' : self._executive_control_echo,
            'F O R C E S   I N   S C A L A R   S P R I N G S        ( C E L A S 4 )' : self._executive_control_echo,

            'F O R C E S   I N   G A P   E L E M E N T S   ( C G A P )' : self._executive_control_echo,

            'F O R C E   D I S T R I B U T I O N   I N   B A R   E L E M E N T S          ( C B A R )' : self._executive_control_echo,
            'F O R C E S   I N   B A R   E L E M E N T S         ( C B A R )' : self._executive_control_echo,
            'F O R C E S   I N   B E A M   E L E M E N T S        ( C B E A M )' : self._executive_control_echo,
            'F O R C E S   I N   B E N D   E L E M E N T S        ( C B E N D )' : self._executive_control_echo,
            'F O R C E S   I N   B U S H   E L E M E N T S        ( C B U S H )' : self._executive_control_echo,

            'F O R C E S   I N   R O D   E L E M E N T S     ( C R O D )' : self._forces_in_crod_elements,
            'F O R C E S   I N   R O D   E L E M E N T S     ( C T U B E )' : self._forces_in_ctube_elements,
            'F O R C E S   I N   R O D   E L E M E N T S     ( C O N R O D )' : self._forces_in_conrod_elements,


            'F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)' : self._executive_control_echo,
            'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = BILIN' : self._forces_in_cquad4s_bilinear,
            'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = CUBIC' : self._executive_control_echo,
            'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
            'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
            'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = CENTER' : self._executive_control_echo,
            'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )        OPTION = SGAGE' : self._executive_control_echo,
            'F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' :  self._executive_control_echo,
            'F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
            'F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,

            # thermal forces
            'H E A T   F L O W   I N   R O D   E L E M E N T S      ( C R O D )' : self._executive_control_echo,
            'H E A T   F L O W   I N   R O D   E L E M E N T S      ( C O N R O D )' : self._executive_control_echo,
            'H E A T   F L O W   I N   R O D   E L E M E N T S      ( C T U B E )' : self._executive_control_echo,

            'H E A T   F L O W   I N   B A R   E L E M E N T S          ( C B A R )' : self._executive_control_echo,
            'H E A T   F L O W   I N   B E A M   E L E M E N T S        ( C B E A M )' : self._executive_control_echo,
            'H E A T   F L O W   I N   B E N D   E L E M E N T S        ( C B E N D )' : self._executive_control_echo,

            'H E A T   F L O W   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            'H E A T   F L O W   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
            'H E A T   F L O W   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
            'H E A T   F L O W   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
            'H E A T   F L O W   I N   A X I S Y M M E T R I C   T R I A N G U L A R   E L E M E N T S   ( C T R I A X 6 )' : self._executive_control_echo,

            'H E A T   F L O W   I N    T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )' : self._executive_control_echo,
            'H E A T   F L O W   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
            'H E A T   F L O W   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._executive_control_echo,

            'H E A T   F L O W   I N T O   H B D Y   E L E M E N T S   (CHBDY)' : self._executive_control_echo,

            # b-list outputs
            'S T R E N G T H   R A T I O S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            'S T R E N G T H   R A T I O S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
            ## ===========================================
            # complex

            # complex stress - not implemented
            'C O M P L E X   S T R E S S E S   I N   S C A L A R   S P R I N G S   ( C E L A S 1 )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   S C A L A R   S P R I N G S   ( C E L A S 2 )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   S C A L A R   S P R I N G S   ( C E L A S 3 )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   S C A L A R   S P R I N G S   ( C E L A S 4 )' : self._executive_control_echo,

            'C O M P L E X   S T R E S S E S   I N   B U S H   E L E M E N T S   ( C B U S H )' : self._executive_control_echo,

            'C O M P L E X   S T R E S S E S   I N   R O D   E L E M E N T S   ( C R O D )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   R O D   E L E M E N T S   ( C T U B E )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   R O D   E L E M E N T S   ( C O N R O D )' : self._executive_control_echo,

            'C O M P L E X   S T R E S S E S   I N   B A R   E L E M E N T S   ( C B A R )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   B E A M   E L E M E N T S   ( C B E A M )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   B E N D   E L E M E N T S   ( C B E N D )' : self._executive_control_echo,

            'C O M P L E X   S T R E S S E S   I N   S H E A R   P A N E L S   ( C S H E A R )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = CENTER' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,

            'C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   H E X A H E D R O N   E L E M E N T S   ( H E X A )' : self._executive_control_echo,

            # complex strain - not implemented
            'C O M P L E X    S T R A I N S    I N   S C A L A R   S P R I N G S   ( C E L A S 1 )' : self._executive_control_echo,
            'C O M P L E X    S T R A I N S    I N   S C A L A R   S P R I N G S   ( C E L A S 2 )' : self._executive_control_echo,
            'C O M P L E X    S T R A I N S    I N   S C A L A R   S P R I N G S   ( C E L A S 3 )' : self._executive_control_echo,

            'C O M P L E X     S T R A I N S   I N   B U S H   E L E M E N T S   ( C B U S H )' : self._executive_control_echo,

            'C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C R O D )' : self._executive_control_echo,
            'C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C T U B E )' : self._executive_control_echo,
            'C O M P L E X    S T R A I N S    I N   R O D   E L E M E N T S   ( C O N R O D )' : self._executive_control_echo,

            'C O M P L E X    S T R A I N S    I N   B A R   E L E M E N T S   ( C B A R )' : self._executive_control_echo,
            'C O M P L E X    S T R A I N S    I N   B E A M   E L E M E N T S   ( C B E A M )' : self._executive_control_echo,
            'C O M P L E X    S T R A I N S    I N   B E N D   E L E M E N T S   ( C B E N D )' : self._executive_control_echo,


            'C O M P L E X     S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            'C O M P L E X     S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
            'C O M P L E X     S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
            'C O M P L E X     S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = CENTER' : self._executive_control_echo,
            'C O M P L E X     S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,
            'C O M P L E X     S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
            'C O M P L E X     S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,

            'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
            'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,

            'C O M P L E X     S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
            'C O M P L E X     S T R A I N S   I N   H E X A H E D R O N   E L E M E N T S   ( H E X A )' : self._executive_control_echo,

            # complex force - not implemented
            'C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 1 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 2 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 3 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   S C A L A R   S P R I N G S   ( C E L A S 4 )' : self._executive_control_echo,

            'C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 1 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 2 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 3 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   S C A L A R   D A M P E R S   ( C D A M P 4 )' : self._executive_control_echo,

            'C O M P L E X   F O R C E S   I N   B U S H   E L E M E N T S   ( C B U S H )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   V I S C   E L E M E N T S   ( C V I S C )' : self._executive_control_echo,

            'C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C R O D )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C T U B E )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   R O D   E L E M E N T S   ( C O N R O D )' : self._executive_control_echo,

            'C O M P L E X   F O R C E S   I N   B A R   E L E M E N T S   ( C B A R )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   B E A M   E L E M E N T S   ( C B E A M )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   B E N D   E L E M E N T S   ( C B E N D )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   B E N D    E L E M E N T S   ( C B E N D )' : self._executive_control_echo,  # why are there different writers for the same table!!!


            'C O M P L E X   F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )        OPTION = CENTER' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,

            # thermal force - not implemented
            'F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S' : self._executive_control_echo,
            'S T R E S S E S   ( F O R C E S )   I N   G A P   E L E M E N T S      ( C G A P )' : self._executive_control_echo,
            ## ===========================================
            # hyperelastic - not implemented
            'S T R E S S E S   I N   H Y P E R E L A S T I C   T R I A N G L E   E L E M E N T S  ( T R I A F D )' : self._executive_control_echo,

            'S T R E S S E S   I N   H Y P E R E L A S T I C   Q U A D R I L A T E R A L   E L E M E N T S  ( QUAD4FD )' : self._executive_control_echo,
            'S T R E S S E S   I N   H Y P E R E L A S T I C   H E X A H E D R O N   E L E M E N T S  ( HEXA8FD )' : self._executive_control_echo,

            # oug nonlinear - not implemented
            'N O N - L I N E A R - F O R C E   V E C T O R' : self._executive_control_echo,

            # stress nonlinear - not implemented
            'N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C R O D )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C T U B E )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S      ( C O N R O D )' : self._executive_control_echo,

            'N O N L I N E A R   S T R E S S E S   I N   B E A M   E L E M E N T S     ( C B E A M )' : self._executive_control_echo,

            'N O N L I N E A R   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S    ( Q U A D 4 )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S    ( Q U A D R )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S      ( T R I A 3 )' : self._executive_control_echo,

            'N O N L I N E A R   S T R E S S E S   I N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S     ( H E X A )' : self._executive_control_echo,



            # nonlinear hyperelastic stress - not implemented
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   T R I A N G L E   E L E M E N T S  ( TRIA3FD )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   Q U A D R I L A T E R A L   E L E M E N T S  ( QUADFD )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   Q U A D R I L A T E R A L   E L E M E N T S  ( QUAD4FD )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S  IN  H Y P E R E L A S T I C   A X I S Y M M.  Q U A D R I L A T E R A L  ELEMENTS (QUADXFD)' : self._executive_control_echo,

            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   T R I A N G L E   E L E M E N T S  ( T R I A F D )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   A X I S Y M M.  T R I A N G L E   E L E M E N T S  (TRIAXFD)' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   A X I S Y M M.  T R I A N G L E   ELEMENTS  (TRIAX3FD)' : self._executive_control_echo,


            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   T E T R A H E D R O N   E L E M E N T S  ( T E T R A F D )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   P E N T A   E L E M E N T S  ( P E N T A F D )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   H E X A H E D R O N   E L E M E N T S  ( H E X A F D )' : self._executive_control_echo,

            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   T E T R A H E D R O N   E L E M E N T S  ( TETRA4FD )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   P E N T A   E L E M E N T S  ( PENTA6FD )' : self._executive_control_echo,
            'N O N L I N E A R   S T R E S S E S   I N   H Y P E R E L A S T I C   H E X A H E D R O N   E L E M E N T S  ( HEXA8FD )' : self._executive_control_echo,

            # nonlinear forces
            'N O N L I N E A R   F O R C E S  A N D  S T R E S S E S  I N   B U S H   E L E M E N T S    ( C B U S H )' : self._executive_control_echo,

            ## ===========================================
            # complex oug - not implemented
            'C O M P L E X   V E L O C I T Y   V E C T O R' : self._executive_control_echo,
            'C O M P L E X   A C C E L E R A T I O N   V E C T O R' : self._executive_control_echo,
            'C O M P L E X   L O A D   V E C T O R' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   O F   S I N G L E   P O I N T   C O N S T R A I N T' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   O F   M U L T I P O I N T   C O N S T R A I N T' : self._executive_control_echo,

            # solution set - not implemented
            'D I S P L A C E M E N T   V E C T O R   (SOLUTION SET)' : self._executive_control_echo,

            # complex solution set - not implemented
            'C O M P L E X   D I S P L A C E M E N T   V E C T O R  (SOLUTION SET)' : self._executive_control_echo,

            # complex weirdness - not implemented
            'C O M P L E X   V E L O C I T I E S   I N   S L O T   E L E M E N T S   ( C S L O T 3 - S T R E S S )' : self._executive_control_echo,
            'C O M P L E X   V E L O C I T I E S   I N   S L O T   E L E M E N T S   ( C S L O T 4 - S T R E S S )' : self._executive_control_echo,
            'C O M P L E X   V E L O C I T I E S   I N   A X I S Y M M E T R I C   F L U I D   E L E M E N T S   ( C A X I F 2 - S T R E S S )' : self._executive_control_echo,
            'C O M P L E X   V E L O C I T I E S   I N   A X I S Y M M E T R I C   F L U I D   E L E M E N T S   ( C A X I F 3 - S T R E S S )' : self._executive_control_echo,
            'C O M P L E X   A C O U S T I C   P R E S S U R E   R E S U L T S' : self._executive_control_echo,

            ## ===========================================
            # p-element stress - not implemented
            'S T R E S S E S   I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )' : self._executive_control_echo,

            'S T R E S S E S   I N   P - V E R S I O N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            'S T R E S S E S   I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,

            'S T R E S S E S   I N   P - V E R S I O N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )' : self._executive_control_echo,
            'S T R E S S E S   I N   P - V E R S I O N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
            'S T R E S S E S   I N   P - V E R S I O N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._executive_control_echo,

            # p-element strain - not implemented
            'S T R A I N S    I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )' : self._executive_control_echo,
            'S T R A I N S    I N   P - V E R S I O N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            'S T R A I N S    I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,

            'S T R A I N S    I N   P - V E R S I O N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._executive_control_echo,
            'S T R A I N S    I N   P - V E R S I O N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
            'S T R A I N S    I N   P - V E R S I O N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )' : self._executive_control_echo,

            # p-element force - not implemented
            'F O R C E S   I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )' : self._executive_control_echo,
            'F O R C E S   I N   P - V E R S I O N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            'F O R C E S   I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,

            # p-element temperature
            'T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S   I N   B E A M   P - E L E M E N T S' : self._executive_control_echo,

            'T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S   I N   T R I A N G U L A R   P - E L E M E N T S' : self._executive_control_echo,
            'T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S   I N   Q U A D R I L A T E R A L   P - E L E M E N T S' : self._executive_control_echo,

            'T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S   I N   T E T R A H E D R O N   P - E L E M E N T S' : self._executive_control_echo,
            'T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S   I N   P E N T A H E D R O N   P - E L E M E N T S' : self._executive_control_echo,
            'T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S   I N   H E X A H E D R O N   P - E L E M E N T S' : self._executive_control_echo,

            ## ===========================================
            # complex p-element stress - not implemented
            'C O M P L E X   S T R E S S E S   I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )' : self._executive_control_echo,

            'C O M P L E X   S T R E S S E S   I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   P - V E R S I O N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,

            'C O M P L E X   S T R E S S E S   I N   P - V E R S I O N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   P - V E R S I O N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
            'C O M P L E X   S T R E S S E S   I N   P - V E R S I O N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )' : self._executive_control_echo,

            # complex p-element strain - not implemented
            'C O M P L E X    S T R A I N S    I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )' : self._executive_control_echo,

            'C O M P L E X    S T R A I N S    I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
            'C O M P L E X    S T R A I N S    I N   P - V E R S I O N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,

            'C O M P L E X    S T R A I N S    I N   P - V E R S I O N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
            'C O M P L E X    S T R A I N S    I N   P - V E R S I O N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )' : self._executive_control_echo,
            'C O M P L E X    S T R A I N S    I N   P - V E R S I O N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._executive_control_echo,

            # complex p-element force - not implemented
            'C O M P L E X   F O R C E S   I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
            'C O M P L E X   F O R C E S   I N   P - V E R S I O N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            ## ===========================================
            # ???
            # peak/rms not included
            'A C C E L E R A T I O N S   V E L O C I T I E S   A N D   P R E S S U R E   L E V E L S' : self._executive_control_echo,

            'F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
            'F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
            'F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,
            'F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
            'F A I L U R E   I N D I C E S   F O R   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,


            # more weirdness - not implemented
            'C O M P L E X   F I E L D S   I N   A A B S F   E L E M E N T S   ( C A A B S F )' : self._executive_control_echo,
            'S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       1' : self._executive_control_echo,
            'S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       2' : self._executive_control_echo,
            'S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       3' : self._executive_control_echo,
            'S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       4' : self._executive_control_echo,
            'S T R E S S E S   A T   G R I D   P O I N T S   - -     S U R F A C E       5' : self._executive_control_echo,
            'R E S U L T S   F O R   S L I D E   L I N E   E L E M E N T S   (IN ELEMENT SYSTEM)' : self._executive_control_echo,
            'P - E L E M E N T    E R R O R    E S T I M A T E    T A B L E' : self._executive_control_echo,
            'P - E L E M E N T    E R R O R    E S T I M A T E    S U M M A R Y    T A B L E' : self._executive_control_echo,
            'S U P E R E L E M E N T   T R E E' : self._executive_control_echo,
            'D E S I G N   C Y C L E =       1    S U B C A S E =       1' : self._executive_control_echo,
            'E L E M E N T   I N T E R N A L   F O R C E S   A N D   M O M E N T S' : self._executive_control_echo,
            'E L E M E N T   S T R E S S   D I S C O N T I N U I T I E S  - -     S U R F A C E       1' : self._executive_control_echo,
            'E L E M E N T   S T R E S S   D I S C O N T I N U I T I E S  - -     S U R F A C E      91' : self._executive_control_echo,
            'P R I N C I P A L   E L E M E N T   S T R E S S   D I S C O N T I N U I T I E S  - -       V O L U M E       92' : self._executive_control_echo,
            'U S E T   D E F I N I T I O N   T A B L E   ( I N T E R N A L   S E Q U E N C E ,   R O W   S O R T )' : self._executive_control_echo,
            'S U M M A T I O N   O F   E L E M E N T   O R I E N T E D   F O R C E S   O N   A D J A C E N T   E L E M E N T S' : self._executive_control_echo,
            '*                  E L E M E N T   P R O P E R T Y   S U M M A R Y     (BY ELEMENT TYPE / ID)              *' : self._executive_control_echo,
            'N O N - L I N E A R   I T E R A T I O N   M O D U L E   S O L U T I O N   D A T A' : self._executive_control_echo,

            #'* * * END OF JOB * * *': self.end(),
        }
        self.markers = self._marker_map.keys()

    def _start_log(self, log=None, debug=False):
        """
        Sets up a dummy logger if one is not provided

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
        i = 0
        while 'PAGE' not in line and '* * * END OF JOB * * *' not in line:
            line = self.infile.readline()[1:].strip()
            lines.append(line)
            self.i += 1
            self.fatal_check(line)
            if not line:
                i += 1
            else:
                i = 0
            if i >= 1000:  # if there are 1000 blank lines, stop
                #print('i=%s line=%r' % (i, line))
                if is_release:
                    break
                else:
                    raise RuntimeError('infinite loop')
        #self.grid_point_weight.read_grid_point_weight(lines)

    def _executive_control_echo(self):
        line = ''
        lines = []
        i = 0
        while 'PAGE' not in line and '* * * END OF JOB * * *' not in line:
            line = self.infile.readline()[1:].strip()
            lines.append(line)
            self.i += 1
            self.fatal_check(line)
            if not line:
                i += 1
            else:
                i = 0
            if i >= 1000:  # if there are 1000 blank lines, stop
                #print('i=%s line=%r' % (i, line))
                if is_release:
                    break
                else:
                    raise RuntimeError('infinite loop')

        #self.grid_point_weight.read_grid_point_weight(lines)

    def fatal_check(self, line):
        if 'FATAL' in line:
            raise FatalError(line)

    def _nastran_file_and_system_parameter_echo(self):
        line = ''
        lines = []
        while 'PAGE' not in line and '* * * END OF JOB * * *' not in line:
            line = self.infile.readline()[1:].strip()
            lines.append(line)
            self.i += 1
            self.fatal_check(line)
        #self.grid_point_weight.read_grid_point_weight(lines)

    def _grid_point_singularity_table(self):
        line = ''
        lines = []
        while 'PAGE' not in line and '* * * END OF JOB * * *' not in line:
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

    def _get_minus_lines(self, debug=False):
        """
        Reads the following section:

        1    MSC.NASTRAN JOB                                                       FEBRUARY  26, 2014  MSC.NASTRAN  6/17/05   PAGE    14
             DEFAULT
        0                                                                                                            SUBCASE 1

                                                             L O A D   V E C T O R

        and finds the lines whose first characters are 1 and 0.  The
        LOAD VECTOR line is line 0 and lines count down from that in self.stored_lines.
        """
        if debug:
            print('-------------------------')
            print(self.stored_lines[0])
        msg = ''
        i0 = None
        i1 = None
        nlines = len(self.stored_lines)
        for i, line in enumerate(self.stored_lines):
            #if line[0].strip():
            if debug:
                print('%-3i %-2i %r' % (i, i - nlines, line[0]))
            if line[0] == '0':
                i0 = i - nlines
            elif line[0] == '1':
                i1 = i - nlines

        # A header block has a 1 for the first character on the line
        # It's closed by a 0.
        # If i1-i0 == -2, -3 (we're working with negative line numbers
        # (from the end of the stack), we can parse it.
        # If we didn't find i1, it means we probably found it before (and
        # since it's repeated, we don't need it) and assume it's a -2 stack.
        # If it's -3, we know this is a nonlinear run.  Our assumption below
        # (assume -2) is invalid for a nonlinear BDF, but that isn't supported
        # anyways.
        #
        # i1 doesn't always exist, but presumably we've found it before,
        # so we don't need to change it
        found_i1 = True
        if i1 is None:
            i1 = i0 - 2
            found_i1 = False
        assert i0 is not None and i1 is not None, 'i1=%s i0=%s'  % (i1, i0)

        # i1 comes before (it's more negative) than i0
        delta = i1 - i0
        if delta not in [-2, -3]:
            if debug:
                msg = 'found_i1=%s i1=%s i0=%s delta=%s nlines=%i; delta should be in [-2, -3]' % (found_i1, i1, i0, delta, nlines)
                raise RuntimeError(msg)
            self.get_minus_lines(debug=True)


        #header_lines = self.stored_lines[i1 : i0+1]
        #assert header_lines[-1][0] == '0'
        #i = -len(header_lines)
        #for i, line in enumerate(header_lines):
            #msg += 'i=%s - %s\n' % (i, line.rstrip()[:118])
        #print(msg)


        header_lines = self.stored_lines[i1 : i0 + 1]
        if debug:
            for i, line in enumerate(header_lines):
                print('  header i=%s*%r' % (i, line.strip()))
        return found_i1, header_lines, delta

    def _read_f06_subcase_header(self, n=-2):
        """
        -4 -> 1                                                     JANUARY   5, 2014  MSC.NASTRAN 11/25/11   PAGE    14
        -3 -> DEFAULT
        -2 -> xxx             subcase 1
        """
        found_i1, header_lines, delta = self._get_minus_lines()
        if found_i1:
            title_line = header_lines[0]
            assert 'PAGE' in title_line, title_line
            self.title = title_line[1:75].strip()
            date = title_line[75:93].strip()
            if date:
                if ',' in date:
                    month, day, year = date.split()
                    assert day[-1] == ',', day
                    day = day[:-1]
                else:
                    raise RuntimeError('Couldnt parse date; line=%r' % title_line.strip())

                self._set_f06_date(month, day, year)
        #self.title = subcaseName  # 'no title'

        subcase_name = ''
        #print("subcase_name = %r" % subcase_name)

        if found_i1:
            subtitle = header_lines[1].strip()
            label_isubcase = header_lines[2]
            assert delta in [-2, -3], delta
        else:
            subtitle = header_lines[1].strip()
            label_isubcase = header_lines[2]
            assert delta == -2, delta
        #print('label_isubcase  = %r' % label_isubcase[:120].strip())
        label, isubcase = self._parse_label_isubcase(label_isubcase, stop_on_failure=False)
        #print('label  = %r' % label.strip())
        #print('isubcase  = %s' % isubcase)

        # TODO: this hardcodes the subtitle, but with bad splitting
        #subtitle = 'subtitle'
        if 'EIGENVALUE' in subtitle:
            msg = 'subtitle=%r has EIGENVALUE in it...\n' % subtitle
            msg += 'subcase = %s\n' % isubcase
            msg += 'label = %s\n' % label
            msg += '\nStored Lines\n'
            msg += '------------\n'

            stored_lines_to_print = self.stored_lines[-10:]
            i = -len(stored_lines_to_print)
            for i, line in enumerate(stored_lines_to_print):
                msg += 'i=%s - %s\n' % (i-10, line.rstrip())
            raise RuntimeError(msg)

        #subtitle = 'SUBCASE %s' % isubcase
        #label = 'SUBCASE %s' % isubcase
        #self.iSubcaseNameMap[self.isubcase] = [self.subtitle, self.label]

#title      date_stamp  page_stamp
#subtitle
#label      ???

        self._subtitle = subtitle
        transient = self.stored_lines[-1].strip()
        is_sort1 = True
        if transient:
            try:
                trans_word, trans_value = transient.split('=')
            except ValueError:
                msg = 'transient = %r\n' % transient
                msg += 'stored lines = [%s]' % ''.join(self.stored_lines)
                #print(msg)
                raise ValueError(msg)
            trans_word = trans_word.strip()
            if ',' in trans_value:
                trans_value = [float(trans_valuei) for trans_valuei in trans_value.split(',')]
            else:
                trans_value = float(trans_value)
            transient = [trans_word, trans_value]

            if trans_word == 'LOAD STEP':  # nonlinear statics
                analysis_code = 10
                is_sort1 = True
            elif trans_word in ['TIME', 'TIME STEP']:  # TODO check name
                analysis_code = 6
                is_sort1 = True
            elif trans_word == 'EIGENVALUE':  # normal modes
                analysis_code = 2
                is_sort1 = True
            elif trans_word == 'FREQ':  # TODO check name
                analysis_code = 5
                is_sort1 = True
            elif trans_word == 'FREQUENCY':
                analysis_code = 5
                is_sort1 = True
            elif trans_word == 'POINT-ID':
                is_sort1 = False
                analysis_code = None
            elif trans_word == 'ELEMENT-ID':
                is_sort1 = False
                analysis_code = None
            elif trans_word == 'COMPLEX EIGENVALUE':
                is_sort1 = True
                analysis_code = 9
            else:
                raise NotImplementedError('transient_word=%r is not supported...' % trans_word)
        else:
            transient = None
            analysis_code = 1


        key = (isubcase, analysis_code, subtitle)
        if key not in self.labels:
            self.subtitles[isubcase].append(subtitle)
            self.labels[key] = label
        self.iSubcaseNameMap[isubcase] = [subtitle, analysis_code, label]

        dt = None
        if transient is not None:
            dt = transient[1]

        if not is_sort1:
            msg = 'trans_word=%r\n' % trans_word
            msg += self.infile.readline().rstrip()
            msg += self.infile.readline().rstrip()
            #msg += self.stored_lines[-1].strip()
            raise RuntimeError(msg)
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
        #eType = eTypeLine[30:40]
        #totalEnergy1 = eTypeLine[99:114]

        mode_line = self.skip(1)[1:]
        iMode = mode_line[24:40]
        #totalEnergy2 = mode_line[99:114]
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

        if not is_release:
            if isubcase in self.iSubcases:
                self.strain_energy[isubcase].readF06Data(data, transient)
            else:
                sed = RealStrainEnergy(data, transient)
                sed.readF06Data(data, transient)
                self.strain_energy[isubcase] = sed

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

            data_code = {
                'analysis_code': analysis_code,
                'device_code': 1, 'sort_code': 0,
                'sort_bits': [0, 0, 0], 'num_wide': 9,

                'table_code': 1,
                'table_name': 'OES1X',

                #'s_code': s_code,
                #'stress_bits': stress_bits,
                #'element_name': element_name, 'element_type': element_type,

                'format_code': 1,
                'nonlinear_factor': dt,
                'data_names':['lsdvmn'],
                'lsdvmn': 1,
                }
        is_sort1 = True
        if isubcase not in self.grid_point_forces:
            self.grid_point_forces[isubcase] = RealGridPointForces(data_code, is_sort1, isubcase, dt)
        self.grid_point_forces[isubcase].add_f06_data(dt, data)
        self.iSubcases.append(isubcase)

    def _read_f06_table(self, Format, debug=False):
        """
        Reads displacement, spc/mpc forces

        :param Format: list of types [int,str,float,float,float] that maps to sline

        .. seealso:: self.parseLine
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

    def parseLine(self, sline, formats):
        """
        :param sline:   list of strings (split line)
        :param formats: list of types [int,str,float,float,float] that maps to sline
        """
        out = []
        for entry, iformat in zip(sline, formats):
            try:
                entry2 = iformat(entry)
            except:
                #print("sline=%r\n entry=%r format=%s" % (sline, entry, iformat))
                #raise
                return None
            out.append(entry2)
        return out

    def _parse_line_blanks(self, sline, formats):
        """allows blanks"""
        out = []

        for entry, iformat in zip(sline, formats):
            if entry.strip():
                try:
                    entry2 = iformat(entry)
                except:
                    print("sline=%r\n entry=%r format=%s" % (sline, entry, formats))
                    raise
            else:
                entry2 = None
                #print("sline=%r\n entry=%r format=%s" % (sline, entry, iformat))
            out.append(entry2)
        return out

    def read_f06(self, f06_filename=None, vectorized=False):
        """
        Reads the F06 file

        :f06_filename: the file to be parsed (None -> GUI)
        """
        self.is_vectorized = vectorized
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

        #if is_binary(f06_filename):
            #raise IOError('f06_filename=%r is a binary file.' % f06_filename)
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
                if '  ' in marker and '   ' not in marker and '=' not in marker:
                    # markers have multiple words in "stupid" format ('CAT' -> 'C A T'),
                    # so spaced tables 'CAT DOG' -> 'C A T  D O G')
                    print("marker = %r" % marker)
                #print('Title  = %r' % self.subtitle)

            if marker in self.markers:
                blank = 0
                #print("\n1*marker = %r" % marker)
                self._marker_map[marker]()
                if self._stop_after_reading_mass and marker in 'O U T P U T   F R O M   G R I D   P O I N T   W E I G H T   G E N E R A T O R':
                    #print('breaking opgwg')
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
                    #print('breaking blank lines')
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
        self.build_vectorization()

    def _process_f06(self):
        #data = [self.disp,self.SpcForces,self.stress,self.isoStress,self.barStress,self.solidStress,self.temperature]
        data_pack = [self.ctetra_stress, self.cpenta_stress, self.chexa_stress,
                     self.ctetra_strain, self.cpenta_strain, self.chexa_strain,]
        for data_set in data_pack:
            for key, data in iteritems(data_set):
                #print('key=%s class=%s' % (key, type(data)))
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
        for i in range(iskip - 1):
            self.infile.readline()
        self.i += iskip
        return self.infile.readline()

    def _get_next_line(self):
        self.i += 1
        return self.infile.readline()

    def _parse_label_isubcase(self, label_isubcase_line, stop_on_failure=True):
        label = label_isubcase_line[1:65].strip()
        isubcase = label_isubcase_line[65:].strip()
        if isubcase:
            isubcase = int(isubcase.split()[-1])
        else:
            if stop_on_failure:
                msg = 'unknown subcase; label_isubcase_line=%r' % label_isubcase_line
                raise RuntimeError(msg)
            else:
                msg = 'unknown subcase; assuming isubcase=1; label_isubcase_line=%r' % label_isubcase_line
                self.log.error(msg)
            isubcase = 1
        return label, isubcase

if __name__ == '__main__':  # pragma: no cover
    from pyNastran.f06.test.test_f06 import main
    main()

