#pylint: disable=C0301,C0111
from __future__ import print_function
#from six import iteritems
#from six.moves import zip, range
#import os
#from codecs import open

#from pyNastran import is_release
#from pyNastran.utils import print_bad_path
#from pyNastran.utils.log import get_logger

#strainEnergyDensity,TemperatureGradientObject
#from pyNastran.op2.tables.oee_energy.oee_objects import RealStrainEnergy

#from pyNastran.f06.tables.oes import OES
#from pyNastran.f06.tables.oug import OUG
#from pyNastran.f06.tables.oqg import OQG
#from pyNastran.f06.tables.oef import OEF
#from pyNastran.f06.tables.lama import LAMA
#from pyNastran.f06.tables.max_min import MAX_MIN
#from pyNastran.f06.f06_writer import F06Writer
#from pyNastran.op2.tables.ogf_gridPointForces.ogf_Objects import RealGridPointForces

#from pyNastran.utils import is_binary_file
#from pyNastran.f06.errors import FatalError


#class F06(OES, OEF, OUG, OQG, LAMA, MAX_MIN, F06Writer):
class F06():
	def __init__(self, debug=False, log=None):
		"""
		Initializes the F06 object

		:param makeGeom:	reads the BDF tables (default=False)
		:param debug:	   prints data about how the F06 was parsed (default=False)
		:param log:		 a logging object to write debug messages to

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
			'N A S T R A N	F I L E	A N D	S Y S T E M	P A R A M E T E R	E C H O' : self._nastran_file_and_system_parameter_echo,
			'N A S T R A N	E X E C U T I V E	C O N T R O L	E C H O' : self._executive_control_echo,
			'C A S E	C O N T R O L	E C H O' : self._case_control_echo,
			'G R I D   P O I N T   S I N G U L A R I T Y   T A B L E' : self._grid_point_singularity_table,

			# dummy
			'E L E M E N T   G E O M E T R Y   T E S T   R E S U L T S   S U M M A R Y' : self._executive_control_echo,
			'M O D E L   S U M M A R Y' : self._executive_control_echo,
			'M O D E L   S U M M A R Y		  BULK = 0' : self._executive_control_echo,
			'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._grid_point_singularity_table,
			'G R I D   P O I N T   F O R C E   B A L A N C E' : self._grid_point_force_balance,
			'N A S T R A N   S O U R C E   P R O G R A M   C O M P I L A T I O N			 SUBDMAP  =  SESTATIC' : self._executive_control_echo,

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
			'OLOAD	RESULTANT' : self._get_oload_resultant,
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
			#'N O N - D I M E N S I O N A L	H I N G E	M O M E N T	D E R I V A T I V E   C O E F F I C I E N T S' : self._nondimensional_hinge_moment_derivative_coeffs,
			#'A E R O S T A T I C   D A T A   R E C O V E R Y   O U T P U T   T A B L E S' : self._aerostatic_data_recovery_output_tables,
			#'S T R U C T U R A L   M O N I T O R   P O I N T   I N T E G R A T E D   L O A D S' : self._structural_monitor_point_integrated_loads,

			'N O N - D I M E N S I O N A L   S T A B I L I T Y   A N D   C O N T R O L   D E R I V A T I V E   C O E F F I C I E N T S' : self._executive_control_echo,
			'N O N - D I M E N S I O N A L	H I N G E	M O M E N T	D E R I V A T I V E   C O E F F I C I E N T S' :  self._executive_control_echo,
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
			'S T R E S S E S   I N   S C A L A R   S P R I N G S		( C E L A S 1 )' : self._stress_in_celas1_elements,
			'S T R E S S E S   I N   S C A L A R   S P R I N G S		( C E L A S 2 )' : self._stress_in_celas2_elements,
			'S T R E S S E S   I N   S C A L A R   S P R I N G S		( C E L A S 3 )' : self._stress_in_celas3_elements,
			'S T R E S S E S   I N   S C A L A R   S P R I N G S		( C E L A S 4 )' : self._stress_in_celas4_elements,

			'S T R E S S E S   I N   R O D   E L E M E N T S	  ( C R O D )' : self._stress_in_crod_elements,
			'S T R E S S E S   I N   R O D   E L E M E N T S	  ( C T U B E )' : self._stress_in_ctube_elements,
			'S T R E S S E S   I N   R O D   E L E M E N T S	  ( C O N R O D )' : self._stress_in_conrod_elements,
			'S T R E S S E S   I N   B A R   E L E M E N T S		  ( C B A R )' : self._stress_in_cbar_elements,

			'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._stress_in_ctria3_elements,
			'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._stress_in_cquad4_elements,
			'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )		OPTION = BILIN' : self._stress_in_cquad4_bilinear_elements,

			'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )': self._stress_in_composite_ctria3_elements,
			'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )': self._stress_in_composite_cquad4_elements,

			'S T R E S S E S   I N	T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )' : self._stress_in_ctetra_elements,
			'S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._stress_in_chexa_elements,
			'S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._stress_in_cpenta_elements,

			# strain - good
			'S T R A I N S	I N   S C A L A R   S P R I N G S		( C E L A S 1 )' : self._strain_in_celas1_elements,
			'S T R A I N S	I N   S C A L A R   S P R I N G S		( C E L A S 2 )' : self._strain_in_celas2_elements,
			'S T R A I N S	I N   S C A L A R   S P R I N G S		( C E L A S 3 )' : self._strain_in_celas3_elements,
			'S T R A I N S	I N   S C A L A R   S P R I N G S		( C E L A S 4 )' : self._strain_in_celas4_elements,

			'S T R A I N S   I N   R O D   E L E M E N T S	  ( C R O D )' : self._strain_in_crod_elements,
			'S T R A I N S   I N   R O D   E L E M E N T S	  ( C O N R O D )' : self._strain_in_conrod_elements,
			'S T R A I N S	I N   B A R   E L E M E N T S		  ( C B A R )' : self._strain_in_cbar_elements,

			'S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._strain_in_ctria3_elements,
			'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._strain_in_cquad4_elements,
			'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )		OPTION = BILIN' : self._strain_in_cquad4_bilinear_elements,

			'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 3 )' : self._strain_in_composite_ctria3_elements,
			'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 4 )' : self._strain_in_composite_cquad4_elements,

			'S T R A I N S   I N	T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )' : self._strain_in_ctetra_elements,
			'S T R A I N S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._strain_in_chexa_elements,
			'S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._strain_in_cpenta_elements,

			# energy - good
			'E L E M E N T   S T R A I N   E N E R G I E S' : self._element_strain_energies,

			#====================================================================

			# stress - not implemented
			'S T R E S S   D I S T R I B U T I O N   I N   B A R   E L E M E N T S	   ( C B A R )' : self._executive_control_echo,
			'S T R E S S E S   I N   B E N D   E L E M E N T S		( C B E N D )' : self._executive_control_echo,
			'S T R E S S E S   I N   B E A M   E L E M E N T S		( C B E A M )' : self._executive_control_echo,

			'S T R E S S E S   ( F O R C E S )   I N   B U S H 1 D   E L E M E N T S   ( C B U S H 1 D )' : self._executive_control_echo,
			'S T R E S S E S   I N   B U S H   E L E M E N T S		( C B U S H )' : self._executive_control_echo,

			'S T R E S S E S   I N   T R I A X 6   E L E M E N T S' : self._executive_control_echo,
			'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
			'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
			'S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,

			'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
			'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )		OPTION = CENTER' : self._executive_control_echo,
			'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )		OPTION = CUBIC' : self._executive_control_echo,
			'S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )		OPTION = SGAGE' : self._executive_control_echo,

			'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
			'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,
			'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
			'S T R E S S E S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,

			# strain - not implemented
			'S T R A I N S   I N   B U S H   E L E M E N T S		( C B U S H )' : self._executive_control_echo,

			'S T R A I N   D I S T R I B U T I O N   I N   B A R   E L E M E N T S	   ( C B A R )' : self._executive_control_echo,
			'S T R A I N S	I N   B E A M   E L E M E N T S		( C B E A M )' : self._executive_control_echo,
			'S T R A I N S	I N   B E N D   E L E M E N T S		( C B E N D )' : self._executive_control_echo,

			'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
			'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
			'S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )		OPTION = CENTER' : self._executive_control_echo,
			'S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
			'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,

			# other real - not implemented
			'S T R E S S E S   I N   U S E R   E L E M E N T S (CDUM8)' : self._executive_control_echo,
			'S T R E S S E S   I N   U S E R   E L E M E N T S (CDUM9)' : self._executive_control_echo,
			'E L E M E N T   S T R A I N   E N E R G I E S   ( A V E R A G E )' : self._executive_control_echo,
			'E L E M E N T   S T R E S S   D I S C O N T I N U I T I E S  - -	 S U R F A C E	  99' : self._executive_control_echo,

			# FORCE - not implemented
			'F O R C E S   I N   S C A L A R   S P R I N G S		( C E L A S 1 )' : self._executive_control_echo,
			'F O R C E S   I N   S C A L A R   S P R I N G S		( C E L A S 2 )' : self._executive_control_echo,
			'F O R C E S   I N   S C A L A R   S P R I N G S		( C E L A S 3 )' : self._executive_control_echo,
			'F O R C E S   I N   S C A L A R   S P R I N G S		( C E L A S 4 )' : self._executive_control_echo,

			'F O R C E S   I N   G A P   E L E M E N T S   ( C G A P )' : self._executive_control_echo,

			'F O R C E   D I S T R I B U T I O N   I N   B A R   E L E M E N T S		  ( C B A R )' : self._executive_control_echo,
			'F O R C E S   I N   B A R   E L E M E N T S		 ( C B A R )' : self._executive_control_echo,
			'F O R C E S   I N   B E A M   E L E M E N T S		( C B E A M )' : self._executive_control_echo,
			'F O R C E S   I N   B E N D   E L E M E N T S		( C B E N D )' : self._executive_control_echo,
			'F O R C E S   I N   B U S H   E L E M E N T S		( C B U S H )' : self._executive_control_echo,

			'F O R C E S   I N   R O D   E L E M E N T S	 ( C R O D )' : self._forces_in_crod_elements,
			'F O R C E S   I N   R O D   E L E M E N T S	 ( C T U B E )' : self._forces_in_ctube_elements,
			'F O R C E S   I N   R O D   E L E M E N T S	 ( C O N R O D )' : self._forces_in_conrod_elements,


			'F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)' : self._executive_control_echo,
			'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )		OPTION = BILIN' : self._forces_in_cquad4s_bilinear,
			'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )		OPTION = CUBIC' : self._executive_control_echo,
			'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
			'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
			'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )		OPTION = CENTER' : self._executive_control_echo,
			'F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )		OPTION = SGAGE' : self._executive_control_echo,
			'F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' :  self._executive_control_echo,
			'F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
			'F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,

			# thermal forces
			'H E A T   F L O W   I N   R O D   E L E M E N T S	  ( C R O D )' : self._executive_control_echo,
			'H E A T   F L O W   I N   R O D   E L E M E N T S	  ( C O N R O D )' : self._executive_control_echo,
			'H E A T   F L O W   I N   R O D   E L E M E N T S	  ( C T U B E )' : self._executive_control_echo,

			'H E A T   F L O W   I N   B A R   E L E M E N T S		  ( C B A R )' : self._executive_control_echo,
			'H E A T   F L O W   I N   B E A M   E L E M E N T S		( C B E A M )' : self._executive_control_echo,
			'H E A T   F L O W   I N   B E N D   E L E M E N T S		( C B E N D )' : self._executive_control_echo,

			'H E A T   F L O W   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
			'H E A T   F L O W   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
			'H E A T   F L O W   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
			'H E A T   F L O W   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
			'H E A T   F L O W   I N   A X I S Y M M E T R I C   T R I A N G U L A R   E L E M E N T S   ( C T R I A X 6 )' : self._executive_control_echo,

			'H E A T   F L O W   I N	T E T R A H E D R O N   S O L I D   E L E M E N T S   ( C T E T R A )' : self._executive_control_echo,
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
			'C O M P L E X   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )		OPTION = CENTER' : self._executive_control_echo,
			'C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
			'C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
			'C O M P L E X   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,

			'C O M P L E X   S T R E S S E S   I N   T E T R A H E D R O N   E L E M E N T S   ( C T E T R A )' : self._executive_control_echo,
			'C O M P L E X   S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
			'C O M P L E X   S T R E S S E S   I N   H E X A H E D R O N   E L E M E N T S   ( H E X A )' : self._executive_control_echo,

			# complex strain - not implemented
			'C O M P L E X	S T R A I N S	I N   S C A L A R   S P R I N G S   ( C E L A S 1 )' : self._executive_control_echo,
			'C O M P L E X	S T R A I N S	I N   S C A L A R   S P R I N G S   ( C E L A S 2 )' : self._executive_control_echo,
			'C O M P L E X	S T R A I N S	I N   S C A L A R   S P R I N G S   ( C E L A S 3 )' : self._executive_control_echo,

			'C O M P L E X	 S T R A I N S   I N   B U S H   E L E M E N T S   ( C B U S H )' : self._executive_control_echo,

			'C O M P L E X	S T R A I N S	I N   R O D   E L E M E N T S   ( C R O D )' : self._executive_control_echo,
			'C O M P L E X	S T R A I N S	I N   R O D   E L E M E N T S   ( C T U B E )' : self._executive_control_echo,
			'C O M P L E X	S T R A I N S	I N   R O D   E L E M E N T S   ( C O N R O D )' : self._executive_control_echo,

			'C O M P L E X	S T R A I N S	I N   B A R   E L E M E N T S   ( C B A R )' : self._executive_control_echo,
			'C O M P L E X	S T R A I N S	I N   B E A M   E L E M E N T S   ( C B E A M )' : self._executive_control_echo,
			'C O M P L E X	S T R A I N S	I N   B E N D   E L E M E N T S   ( C B E N D )' : self._executive_control_echo,


			'C O M P L E X	 S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
			'C O M P L E X	 S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
			'C O M P L E X	 S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
			'C O M P L E X	 S T R A I N S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )		OPTION = CENTER' : self._executive_control_echo,
			'C O M P L E X	 S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,
			'C O M P L E X	 S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
			'C O M P L E X	 S T R A I N S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,

			'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
			'S T R A I N S   I N   L A Y E R E D   C O M P O S I T E   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,

			'C O M P L E X	 S T R A I N S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
			'C O M P L E X	 S T R A I N S   I N   H E X A H E D R O N   E L E M E N T S   ( H E X A )' : self._executive_control_echo,

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
			'C O M P L E X   F O R C E S   I N   B E N D	E L E M E N T S   ( C B E N D )' : self._executive_control_echo,  # why are there different writers for the same table!!!


			'C O M P L E X   F O R C E S   A C T I N G   O N   S H E A R   P A N E L   E L E M E N T S   (CSHEAR)' : self._executive_control_echo,
			'C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
			'C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 8 )' : self._executive_control_echo,
			'C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )' : self._executive_control_echo,
			'C O M P L E X   F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D R )		OPTION = CENTER' : self._executive_control_echo,
			'C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
			'C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A 6 )' : self._executive_control_echo,
			'C O M P L E X   F O R C E S   I N   T R I A N G U L A R   E L E M E N T S   ( T R I A R )' : self._executive_control_echo,

			# thermal force - not implemented
			'F I N I T E   E L E M E N T   T E M P E R A T U R E   G R A D I E N T S   A N D   F L U X E S' : self._executive_control_echo,
			'S T R E S S E S   ( F O R C E S )   I N   G A P   E L E M E N T S	  ( C G A P )' : self._executive_control_echo,
			## ===========================================
			# hyperelastic - not implemented
			'S T R E S S E S   I N   H Y P E R E L A S T I C   T R I A N G L E   E L E M E N T S  ( T R I A F D )' : self._executive_control_echo,

			'S T R E S S E S   I N   H Y P E R E L A S T I C   Q U A D R I L A T E R A L   E L E M E N T S  ( QUAD4FD )' : self._executive_control_echo,
			'S T R E S S E S   I N   H Y P E R E L A S T I C   H E X A H E D R O N   E L E M E N T S  ( HEXA8FD )' : self._executive_control_echo,

			# oug nonlinear - not implemented
			'N O N - L I N E A R - F O R C E   V E C T O R' : self._executive_control_echo,

			# stress nonlinear - not implemented
			'N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S	  ( C R O D )' : self._executive_control_echo,
			'N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S	  ( C T U B E )' : self._executive_control_echo,
			'N O N L I N E A R   S T R E S S E S   I N   R O D   E L E M E N T S	  ( C O N R O D )' : self._executive_control_echo,

			'N O N L I N E A R   S T R E S S E S   I N   B E A M   E L E M E N T S	 ( C B E A M )' : self._executive_control_echo,

			'N O N L I N E A R   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S	( Q U A D 4 )' : self._executive_control_echo,
			'N O N L I N E A R   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S	( Q U A D R )' : self._executive_control_echo,
			'N O N L I N E A R   S T R E S S E S   I N   T R I A N G U L A R   E L E M E N T S	  ( T R I A 3 )' : self._executive_control_echo,

			'N O N L I N E A R   S T R E S S E S   I N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )' : self._executive_control_echo,
			'N O N L I N E A R   S T R E S S E S   I N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
			'N O N L I N E A R   S T R E S S E S   I N   H E X A H E D R O N   S O L I D   E L E M E N T S	 ( H E X A )' : self._executive_control_echo,



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
			'N O N L I N E A R   F O R C E S  A N D  S T R E S S E S  I N   B U S H   E L E M E N T S	( C B U S H )' : self._executive_control_echo,

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
			'S T R A I N S	I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )' : self._executive_control_echo,
			'S T R A I N S	I N   P - V E R S I O N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,
			'S T R A I N S	I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,

			'S T R A I N S	I N   P - V E R S I O N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._executive_control_echo,
			'S T R A I N S	I N   P - V E R S I O N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
			'S T R A I N S	I N   P - V E R S I O N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )' : self._executive_control_echo,

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
			'C O M P L E X	S T R A I N S	I N   P - V E R S I O N   B E A M   E L E M E N T S   ( B E A M )' : self._executive_control_echo,

			'C O M P L E X	S T R A I N S	I N   P - V E R S I O N   T R I A N G U L A R   E L E M E N T S   ( T R I A 3 )' : self._executive_control_echo,
			'C O M P L E X	S T R A I N S	I N   P - V E R S I O N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )' : self._executive_control_echo,

			'C O M P L E X	S T R A I N S	I N   P - V E R S I O N   P E N T A H E D R O N   S O L I D   E L E M E N T S   ( P E N T A )' : self._executive_control_echo,
			'C O M P L E X	S T R A I N S	I N   P - V E R S I O N   T E T R A H E D R O N   S O L I D   E L E M E N T S   ( T E T R A )' : self._executive_control_echo,
			'C O M P L E X	S T R A I N S	I N   P - V E R S I O N   H E X A H E D R O N   S O L I D   E L E M E N T S   ( H E X A )' : self._executive_control_echo,

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
			'S T R E S S E S   A T   G R I D   P O I N T S   - -	 S U R F A C E	   1' : self._executive_control_echo,
			'S T R E S S E S   A T   G R I D   P O I N T S   - -	 S U R F A C E	   2' : self._executive_control_echo,
			'S T R E S S E S   A T   G R I D   P O I N T S   - -	 S U R F A C E	   3' : self._executive_control_echo,
			'S T R E S S E S   A T   G R I D   P O I N T S   - -	 S U R F A C E	   4' : self._executive_control_echo,
			'S T R E S S E S   A T   G R I D   P O I N T S   - -	 S U R F A C E	   5' : self._executive_control_echo,
			'R E S U L T S   F O R   S L I D E   L I N E   E L E M E N T S   (IN ELEMENT SYSTEM)' : self._executive_control_echo,
			'P - E L E M E N T	E R R O R	E S T I M A T E	T A B L E' : self._executive_control_echo,
			'P - E L E M E N T	E R R O R	E S T I M A T E	S U M M A R Y	T A B L E' : self._executive_control_echo,
			'S U P E R E L E M E N T   T R E E' : self._executive_control_echo,
			'D E S I G N   C Y C L E =	   1	S U B C A S E =	   1' : self._executive_control_echo,
			'E L E M E N T   I N T E R N A L   F O R C E S   A N D   M O M E N T S' : self._executive_control_echo,
			'E L E M E N T   S T R E S S   D I S C O N T I N U I T I E S  - -	 S U R F A C E	   1' : self._executive_control_echo,
			'E L E M E N T   S T R E S S   D I S C O N T I N U I T I E S  - -	 S U R F A C E	  91' : self._executive_control_echo,
			'P R I N C I P A L   E L E M E N T   S T R E S S   D I S C O N T I N U I T I E S  - -	   V O L U M E	   92' : self._executive_control_echo,
			'U S E T   D E F I N I T I O N   T A B L E   ( I N T E R N A L   S E Q U E N C E ,   R O W   S O R T )' : self._executive_control_echo,
			'S U M M A T I O N   O F   E L E M E N T   O R I E N T E D   F O R C E S   O N   A D J A C E N T   E L E M E N T S' : self._executive_control_echo,
			'*				  E L E M E N T   P R O P E R T Y   S U M M A R Y	 (BY ELEMENT TYPE / ID)			  *' : self._executive_control_echo,
			'N O N - L I N E A R   I T E R A T I O N   M O D U L E   S O L U T I O N   D A T A' : self._executive_control_echo,

			#'* * * END OF JOB * * *': self.end(),
		}
		self.markers = self._marker_map.keys()

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
		  POINT	TYPE   FAILED	  STIFFNESS	   OLD USET		   NEW USET
		   ID			DIRECTION	  RATIO	 EXCLUSIVE  UNION   EXCLUSIVE  UNION
			1		G	  4		 0.00E+00		  B		F		 SB	   S	*
			1		G	  5		 0.00E+00		  B		F		 SB	   S	*
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

		1	MSC.NASTRAN JOB													   FEBRUARY  26, 2014  MSC.NASTRAN  6/17/05   PAGE	14
			 DEFAULT
		0																											SUBCASE 1

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
		-4 -> 1													 JANUARY   5, 2014  MSC.NASTRAN 11/25/11   PAGE	14
		-3 -> DEFAULT
		-2 -> xxx			 subcase 1
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

#title	  date_stamp  page_stamp
#subtitle
#label	  ???

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

				  ELEMENT-TYPE = QUAD4			   * TOTAL ENERGY OF ALL ELEMENTS IN PROBLEM	 =  -1.188367E-05
					 MODE			   1			* TOTAL ENERGY OF ALL ELEMENTS IN SET	  -1 =  -1.188367E-05

									  ELEMENT-ID		  STRAIN-ENERGY		   PERCENT OF TOTAL	STRAIN-ENERGY-DENSITY
											   1		 -5.410134E-08				-0.0929			 -4.328107E-05
											   2		 -3.301516E-09				-0.0057			 -2.641213E-06
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
		self.process()
		if hasattr(self, '_ieigenvalue'):
			del self._ieigenvalue
		self.build_vectorization()

	def process(self):
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

	def make_stamp(Title, today=None):
		if 'Title' is None:
			Title = ''

		#lenghts = [7, 8, 5, 5, 3, 4, 4, 6, 9, 7, 8, 8]
		months = [' January', 'February', 'March', 'April', 'May', 'June',
				  'July', 'August', 'September', 'October', 'November', 'December']
		if today is None:
			today = date.today()
			str_month = months[today.month - 1].upper()
			str_today = '%-9s %2s, %4s' % (str_month, today.day, today.year)
		else:
			(month, day, year) = today
			str_month = months[month - 1].upper()
			str_today = '%-9s %2s, %4s' % (str_month, day, year)
		str_today = str_today  #.strip()

		release_date = '02/08/12'  # pyNastran.__releaseDate__
		release_date = ''
		build = 'pyNastran v%s %s' % (pyNastran.__version__, release_date)
		if Title is None:
			Title = ''
		out = '1	%-67s   %-19s %-22s PAGE %%5i\n' % (Title.strip(), str_today, build)
		return out


	def make_f06_header():
		spaces = ''
		lines1 = [
			spaces + '/* -------------------------------------------------------------------  */\n',
			spaces + '/*							  PYNASTRAN							   */\n',
			spaces + '/*					  - NASTRAN FILE INTERFACE -					  */\n',
			spaces + '/*																	  */\n',
			spaces + '/*			  A Python reader/editor/writer for the various		   */\n',
			spaces + '/*						NASTRAN file formats.						 */\n',
			spaces + '/*					   Copyright (C) 2011-2013						*/\n',
			spaces + '/*			   Steven Doyle, Al Danial, Marcin Garrozik			   */\n',
			spaces + '/*																	  */\n',
			spaces + '/*	This program is free software; you can redistribute it and/or	 */\n',
			spaces + '/*	modify it under the terms of the GNU Lesser General Public		*/\n',
			spaces + '/*	License as published by the Free Software Foundation;			 */\n',
			spaces + '/*	version 3 of the License.										 */\n',
			spaces + '/*																	  */\n',
			spaces + '/*	This program is distributed in the hope that it will be useful,   */\n',
			spaces + '/*	but WITHOUT ANY WARRANTY; without even the implied warranty of	*/\n',
			spaces + '/*	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the	  */\n',
			spaces + '/*	GNU Lesser General Public License for more details.			   */\n',
			spaces + '/*																	  */\n',
			spaces + '/*	You should have received a copy of the GNU Lesser General Public  */\n',
			spaces + '/*	License along with this program; if not, write to the			 */\n',
			spaces + '/*	Free Software Foundation, Inc.,								   */\n',
			spaces + '/*	675 Mass Ave, Cambridge, MA 02139, USA.						   */\n',
			spaces + '/* -------------------------------------------------------------------  */\n',
			'\n']

		spaces = 46 * ' '
		version = 'Version %8s' % pyNastran.__version__
		lines2 = [
			spaces + '* * * * * * * * * * * * * * * * * * * *\n',
			spaces + '* * * * * * * * * * * * * * * * * * * *\n',
			spaces + '* *								 * *\n',
			spaces + '* *								 * *\n',
			spaces + '* *								 * *\n',
			spaces + '* *								 * *\n',
			spaces + '* *			pyNastran			* *\n',
			spaces + '* *								 * *\n',
			spaces + '* *								 * *\n',
			spaces + '* *								 * *\n',
			spaces + '* *%s* *\n' % version.center(33),
			spaces + '* *								 * *\n',
			spaces + '* *								 * *\n',
			spaces + '* *		  %15s		* *\n' % pyNastran.__releaseDate2__,
			spaces + '* *								 * *\n',
			spaces + '* *			Questions			* *\n',
			spaces + '* *		mesheb82@gmail.com	   * *\n',
			spaces + '* *								 * *\n',
			spaces + '* *								 * *\n',
			spaces + '* *								 * *\n',
			spaces + '* * * * * * * * * * * * * * * * * * * *\n',
			spaces + '* * * * * * * * * * * * * * * * * * * *\n\n\n']
		return ''.join(lines1 + lines2)


	def sorted_bulk_data_header():
		msg = '0												 S O R T E D   B U L K   D A T A   E C H O										 \n'
		msg += '				 ENTRY																											  \n'
		msg += '				 COUNT		.   1  ..   2  ..   3  ..   4  ..   5  ..   6  ..   7  ..   8  ..   9  ..  10  .					  \n'
		return msg


	def make_end(end_flag=False):
		lines = []
		lines2 = []
		if end_flag:
			lines = [
				'', '',
				'0								   * * * *  A N A L Y S I S  S U M M A R Y  T A B L E  * * * *',
				'0 SEID  PEID PROJ VERS APRCH	  SEMG SEMR SEKR SELG SELR MODES DYNRED SOLLIN PVALID SOLNL LOOPID DESIGN CYCLE SENSITIVITY',
				' --------------------------------------------------------------------------------------------------------------------------']
			#0	 0	1	1 '		'	T	T	T	T	T	 F	  F	  T	  0	 F	 -1			0		   F

			seid = 0
			peid = 0
			proj = 1
			vers = 1
			approach = '		'

			SELG = 'T'
			SEMG = 'T'
			SEMR = 'T'
			SEKR = 'T'
			SELR = 'T'
			MODES = 'F'
			DYNRED = 'F'

			SOLLIN = 'T'
			PVALID = 0
			SOLNL = 'F'
			LOOPID = -1
			CYCLE = 0
			SENSITIVITY = 'F'

			msg = '	 %s	 %s	%s	%s %8r	%s	%s	%s	%s	%s	 %s	  %s	  %s	  %s	 %s	 %s			%s		   %s' % (
				seid, peid, proj, vers, approach, SEMG, SEMR, SEKR, SELG, SELR, MODES, DYNRED, SOLLIN, PVALID, SOLNL,
				LOOPID, CYCLE, SENSITIVITY)
			lines.append(msg)

			lines2 = [
				'0SEID = SUPERELEMENT ID.',
				' PEID = PRIMARY SUPERELEMENT ID OF IMAGE SUPERELEMENT.',
				' PROJ = PROJECT ID NUMBER.',
				' VERS = VERSION ID.',
				' APRCH = BLANK FOR STRUCTURAL ANALYSIS.  HEAT FOR HEAT TRANSFER ANALYSIS.',
				' SEMG = STIFFNESS AND MASS MATRIX GENERATION STEP.',
				' SEMR = MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES).',
				' SEKR = STIFFNESS MATRIX REDUCTION STEP.',
				' SELG = LOAD MATRIX GENERATION STEP.',
				' SELR = LOAD MATRIX REDUCTION STEP. ',
				' MODES = T (TRUE) IF NORMAL MODES OR BUCKLING MODES CALCULATED.',
				' DYNRED = T (TRUE) MEANS GENERALIZED DYNAMIC AND/OR COMPONENT MODE REDUCTION PERFORMED.',
				' SOLLIN = T (TRUE) IF LINEAR SOLUTION EXISTS IN DATABASE.',
				' PVALID = P-DISTRIBUTION ID OF P-VALUE FOR P-ELEMENTS',
				' LOOPID = THE LAST LOOPID VALUE USED IN THE NONLINEAR ANALYSIS.  USEFUL FOR RESTARTS.',
				' SOLNL = T (TRUE) IF NONLINEAR SOLUTION EXISTS IN DATABASE.',
				' DESIGN CYCLE = THE LAST DESIGN CYCLE (ONLY VALID IN OPTIMIZATION).',
				' SENSITIVITY = SENSITIVITY MATRIX GENERATION FLAG.',
				' ',
				' No PARAM values were set in the Control File.'
			]

		lines3 = [
			' ',
			'1										* * * END OF JOB * * *',
			' ',
			' '
		]
		return '\n'.join(lines + lines2 + lines3)


	#def __init__(self):
		#OP2_F06_Common.__init__(self)
		#self.card_count = {}
		#self.additional_matrices = {}
		#self.matrices = {}
		#self.subcase_key = defaultdict(list)
		#self._results = ResultSet(self.get_all_results())

	def get_all_results(self):
		all_results = ['stress', 'strain', 'element_forces', 'constraint_forces'] + self.get_table_types()
		return all_results

	def _clear_results(self):
		self._results.clear()

	def add_results(self, results):
		if isinstance(results, string_types):
			results = [results]
		all_results = self.get_all_results()
		for result in results:
			result = str(result)
			if result not in all_results:
				raise RuntimeError('%r is not a valid result to remove; all_results=%s' % (result, all_results))
			if 'stress' == result:
				stress_results = []
				for result in all_results:
					if 'stress' in result.lower():
						stress_results.append(result)
				#stress_results = [result if 'stress' in result.lower() for result in all_results]
				self._results.update(stress_results)
			elif 'strain' == result:
				strain_results = []
				for result in all_results:
					if 'strain' in result.lower():
						strain_results.append(result)
				#strain_results = [result if 'strain' in result.lower() for result in all_results]
				self._results.update(strain_results)
			elif 'stress' in result.lower():
				self._results.add('stress')
			elif 'strain' in result.lower():
				self._results.add('strain')
			elif result in ('spc_forces', 'mpc_forces', 'constraint_forces'):
				self._results.add('constraint_forces')
			elif 'force' in result.lower(): # could use more validation...
				self._results.add('element_forces')
			# thermalLoad_VU_3D, thermalLoad_1D, thermalLoad_CONV, thermalLoad_2D_3D
			self._results.add(result)

	def set_results(self, results):
		if isinstance(results, string_types):
			results = [results]
		self._clear_results()
		self.add_results(results)

	def remove_results(self, results):
		self._results.remove(results)

	def make_f06_header(self):
		"""If this class is inherited, the F06 Header may be overwritten"""
		return make_f06_header()

	def make_stamp(self, Title, today):
		"""If this class is inherited, the PAGE stamp may be overwritten"""
		return make_stamp(Title, today)

	def make_grid_point_singularity_table(self, failed):
		msg = ''
		if failed:
			msg += '0										 G R I D   P O I N T   S I N G U L A R I T Y   T A B L E\n'
			msg += '0							 POINT	TYPE   FAILED	  STIFFNESS	   OLD USET		   NEW USET\n'
			msg += '							   ID			DIRECTION	  RATIO	 EXCLUSIVE  UNION   EXCLUSIVE  UNION\n'
			for (nid, dof) in failed:
				msg += '						 %8s		G	  %s		 0.00E+00		  B		F		 SB	   SB   *\n' % (nid, dof)
		else:
			msg += 'No constraints have been applied...\n'

		page_stamp = self.make_stamp(self.title, self.date)
		msg += page_stamp % self.page_num
		self.page_num += 1
		return msg

	def _write_summary(self, f06, card_count=None):
		summary_header = '										M O D E L   S U M M A R Y\n\n'
		summary = ''

		self.cards_to_read = set([

			# rigid elements
			'RBAR', 'RBAR1', 'RBE1', 'RBE2', 'RBE3',

			# spc/mpc constraints
			'SPC', 'SPCADD', 'SPC1', 'SPCD', 'SPCAX',
			'MPC', 'MPCADD',
			'SUPORT', 'SUPORT1',

			# aero cards
			'CAERO1', 'CAERO2', 'CAERO3', 'CAERO4', 'CAERO5',

			# temperature cards
			'CHBDYE', 'CHBDYG', 'CHBDYP',
			'CONV',
		])


		blocks = [
			['POINTS', ['GRID', 'GRDSET', ]],
			['ENTRIES', ['SPOINT']],

			['ELEMENTS',
			 [
				 # these are sorted
				 # elements
				 'CONM1', 'CONM2', 'CMASS1', 'CMASS2', 'CMASS3', 'CMASS4',

				 # springs
				 'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4', 'CELAS5',

				 # bushings
				 'CBUSH', 'CBUSH1D', 'CBUSH2D',

				 # dampers
				 'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4', 'CDAMP5',

				 # bar flags
				 'BAROR', 'CBARAO',
				 # bars
				 'CBAR', 'CROD', 'CTUBE', 'CBEAM', 'CBEAM3', 'CONROD', 'CBEND',

				 # shells
				 'CTRIA3', 'CTRIA6', 'CTRIAR', 'CTRIAX', 'CTRIAX6',
				 'CQUAD4', 'CQUAD8', 'CQUADR', 'CQUADX', 'CQUAD',

				 # solids
				 'CTETRA', 'CPENTA', 'CHEXA',

				 # other
				 'CSHEAR', 'CVISC', 'CRAC2D', 'CRAC3D',
				 'CGAP', 'CFAST', 'RBE2', 'RBE3',

				 # thermal
				 'CHBDYP', 'CHBDYG', 'CONV',
			 ]],
		]
		#print("self.card_count", self.card_count)
		if card_count is None:
			card_count = self.card_count

		for block in blocks:
			block_name, keys = block
			key_count = 0
			for key in sorted(keys):
				try:
					value = card_count[key]
					summary += '								   NUMBER OF %-8s %-8s = %8s\n' % (key, block_name, value)
					key_count += 1
				except KeyError:
					pass
			if key_count:
				summary += ' \n'

		if summary:
			f06.write(summary_header)
			f06.write(summary)

			page_stamp = self.make_stamp(self.title, self.date)
			f06.write(page_stamp % self.page_num)
			self.page_num += 1

	def write_f06(self, f06_outname, is_mag_phase=False, is_sort1=True,
				  delete_objects=True, end_flag=False, quiet=False):
		"""
		Writes an F06 file based on the data we have stored in the object

		Parameters
		----------
		f06_outname : str
			the name of the F06 file to write
		is_mag_phase : bool; default=False
			should complex data be written using Magnitude/Phase
			instead of Real/Imaginary
			Real objects don't use this parameter
		is_sort1 : bool; default=True
			writes output in SORT1 format if the output is transient;
			ignored for static analyses
		delete_objects : bool; default=True
			should objects be deleted after they're written to reduce memory
		end_flag : bool; default=False
			should a dummy Nastran "END" table be made
		quiet : bool; default=False
			 suppress print messages
		"""
		if not quiet:
			print("F06:")
		if isinstance(f06_outname, str):
			if PY2:
				f06 = open(f06_outname, 'wb')
			else:
				f06 = open(f06_outname, 'w')
			self._write_summary(f06)
		else:
			assert isinstance(f06_outname, file), 'type(f06_outname)= %s' % f06_outname
			f06 = f06_outname
			f06_outname = f06.name
			print('f06_outname =', f06_outname)

		page_stamp = self.make_stamp(self.title, self.date)
		if self.grid_point_weight.reference_point is not None:
			if not quiet:
				print("grid_point_weight")
			self.page_num = self.grid_point_weight.write_f06(f06, page_stamp, self.page_num)
			assert isinstance(self.page_num, int), self.grid_point_weight.__class__.__name__

		if self.oload_resultant is not None:
			self.page_num = self.oload_resultant.write_f06(f06, page_stamp, self.page_num)
			assert isinstance(self.page_num, int), self.oload_resultant.__class__.__name__

		# writes all results for
		self._write_f06_subcase_based(f06, page_stamp, delete_objects=delete_objects,
									  is_mag_phase=is_mag_phase, is_sort1=is_sort1,
									  quiet=quiet)
		#self._write_f06_time_based(f06, page_stamp)
		f06.write(make_end(end_flag))
		f06.close()

	def _write_f06_subcase_based(self, f06, page_stamp, delete_objects=True,
								 is_mag_phase=False, is_sort1=True, quiet=False):
		header = ['	 DEFAULT																														\n',
				  '\n', '']

		# eigenvalues are written first
		f06.write(page_stamp % self.page_num)
		self.page_num += 1
		for ikey, result in sorted(iteritems(self.eigenvalues)):
			if not quiet:
				print('%-18s case=%r' % (result.__class__.__name__, ikey))
			self.page_num = result.write_f06(f06, header, page_stamp,
											 page_num=self.page_num)
			assert isinstance(self.page_num, int), 'pageNum=%r' % str(self.page_num)
			if delete_objects:
				del result
			self.page_num += 1

		# then eigenvectors
		# has a special header
		# isubcases = sorted(self.iSubcaseNameMap.keys())

		# TODO: superelement version...need the nominal...
		res_keys_subcase = self.subcase_key

		for isubcase, res_keys in sorted(iteritems(res_keys_subcase)):
			for res_key in res_keys:
				if isinstance(res_key, tuple):
					is_compressed = False
				else:
					# int
					is_compressed = True
					isubcase = res_key

				if res_key not in self.eigenvectors:
					continue
				result = self.eigenvectors[res_key]
				subtitle = result.subtitle
				header[0] = '	 %s\n' % subtitle
				header[1] = '0																											SUBCASE %i\n' % isubcase
				#header[2] = complex/nonlinear

				res_length = 18
				res_format = '*%%-%is SUBCASE=%%i' % res_length
				res_format_vectorized = ' %%-%is SUBCASE=%%i SUBTITLE=%%s' % res_length
				class_name = result.__class__.__name__
				if hasattr(result, 'data'):
					if not quiet:
						print(res_format_vectorized % (class_name, isubcase, subtitle))
				else:
					print(res_format % (class_name, isubcase))

				self.page_num = result.write_f06(f06, header, page_stamp,
												 self.page_num, is_mag_phase=is_mag_phase, is_sort1=True)
				assert isinstance(self.page_num, int), 'pageNum=%r' % str(self.page_num)
				if delete_objects:
					del result
				self.page_num += 1

		# finally, we writte all the other tables
		# nastran puts the tables in order of the Case Control deck,
		# but we're lazy so we just hardcode the order

		# subcase name, subcase ID, transient word & value
		header_old = ['	 DEFAULT																														\n',
					  '\n', ' \n']
		header = copy.deepcopy(header_old)
		res_types = [
			self.displacements, self.displacementsPSD, self.displacementsATO, self.displacementsRMS,
			self.displacements_scaled,  # ???
			self.accelerations, self.accelerationsPSD,

			self.force_vectors,
			self.load_vectors,
			self.temperatures,
			self.velocities, self.velocitiesPSD,
			#self.eigenvectors,
			self.eigenvectors_RADCONS,
			self.eigenvectors_RADEFFM,
			self.eigenvectors_RADEATC,

			self.mpc_forces, self.mpc_forcesPSD, self.mpc_forcesATO, self.mpc_forcesRMS,
			self.mpc_forces_RAQCONS,
			#self.mpc_forces_RAQEATC,

			self.spc_forces, self.spc_forcesPSD, self.spc_forcesATO, self.spc_forcesRMS,
			self.thermal_load_vectors,

			#self.strain_energy,
			self.cquad4_strain_energy, self.cquad8_strain_energy,
			self.cquadr_strain_energy, self.cquadx_strain_energy,
			self.ctria3_strain_energy, self.ctria6_strain_energy,
			self.ctriar_strain_energy, self.ctriax_strain_energy,
			self.ctriax6_strain_energy,
			self.ctetra_strain_energy, self.cpenta_strain_energy,
			self.chexa_strain_energy, self.cpyram_strain_energy,
			self.crod_strain_energy, self.ctube_strain_energy,
			self.conrod_strain_energy,
			self.cbar_strain_energy, self.cbeam_strain_energy,
			self.cgap_strain_energy, self.cbush_strain_energy,
			self.celas1_strain_energy, self.celas2_strain_energy,
			self.celas3_strain_energy, self.celas4_strain_energy,
			self.cdum8_strain_energy, self.dmig_strain_energy,
			self.cbend_strain_energy,
			self.genel_strain_energy, self.cshear_strain_energy,
			#------------------------------------------
			# OEF - forces

			# alphabetical order...
			# bars
			self.cbar_force,
			self.cbar_force_10nodes,

			# beam
			self.cbend_force,
			self.cbeam_force,

			# alphabetical
			self.celas1_force,
			self.celas2_force,
			self.celas3_force,
			self.celas4_force,

			self.cquad4_force,
			self.cquad8_force,
			self.cquadr_force,

			self.conrod_force,
			self.crod_force,
			self.cshear_force,
			self.ctria3_force,
			self.ctria6_force,
			self.ctriar_force,
			self.ctube_force,

			# springs
			self.celas1_force,
			self.celas2_force,
			self.celas3_force,
			self.celas4_force,

			# dampers
			self.cdamp1_force,
			self.cdamp2_force,
			self.cdamp3_force,
			self.cdamp4_force,

			# other
			self.cbush_force,
			self.cgap_force,
			self.cvisc_force,

			self.chexa_pressure_force,
			self.cpenta_pressure_force,
			self.ctetra_pressure_force,

			self.coneax_force,

			#------------------------------------------
			# OES - strain
			# 1.  cbar
			# 2.  cbeam
			# 3.  crod/ctube/conrod

			# springs
			self.celas1_strain,
			self.celas2_strain,
			self.celas3_strain,
			self.celas4_strain,

			self.nonlinear_celas1_stress,
			self.nonlinear_celas3_stress,

			# bars/beams
			self.cbar_strain,
			self.cbar_strain_10nodes,
			self.cbeam_strain,

			# plates
			self.cquad4_composite_strain,
			self.cquad8_composite_strain,
			self.cquadr_composite_strain,
			self.ctria3_composite_strain,
			self.ctria6_composite_strain,
			self.ctriar_composite_strain,

			self.nonlinear_ctria3_strain,
			self.nonlinear_cquad4_strain,
			self.ctriax_strain,

			# rods
			self.nonlinear_crod_strain,
			self.nonlinear_ctube_strain,
			self.nonlinear_conrod_strain,

			self.chexa_strain,
			self.conrod_strain,
			self.cpenta_strain,
			self.cquad4_strain,
			self.cquad8_strain,
			self.cquadr_strain,
			self.crod_strain,
			self.cshear_strain,
			self.ctetra_strain,
			self.ctria3_strain,
			self.ctria6_strain,
			self.ctriar_strain,
			self.ctube_strain,

			# bush
			self.cbush_strain,
			self.nonlinear_cbush_stress,
			self.cbush1d_stress_strain,
			#------------------------------------------
			# cbars/cbeams
			self.cbar_stress,
			self.cbar_stress_10nodes,
			self.nonlinear_cbeam_stress,
			self.cbeam_stress,

			# bush
			self.cbush_stress,

			# rods
			self.nonlinear_crod_stress,
			self.nonlinear_ctube_stress,
			self.nonlinear_conrod_stress,

			# shear
			# OES - stress
			self.celas1_stress,
			self.celas2_stress,
			self.celas3_stress,
			self.celas4_stress,

			self.chexa_stress,
			self.conrod_stress,
			self.cpenta_stress,
			self.cquad4_stress,
			self.cquad8_stress,
			self.cquadr_stress,

			self.crod_stress,
			self.cshear_stress,
			self.ctetra_stress,
			self.ctria3_stress,
			self.ctria6_stress,
			self.ctriar_stress,
			self.ctube_stress,

			self.cquad4_composite_stress,
			self.cquad8_composite_stress,
			self.cquadr_composite_stress,
			self.ctria3_composite_stress,
			self.ctria6_composite_stress,
			self.ctriar_composite_stress,

			self.nonlinear_ctria3_stress,
			self.nonlinear_cquad4_stress,
			self.ctriax_stress,

			self.hyperelastic_cquad4_strain,

			#------------------------------------------
			#OEF - Fluxes - tCode=4 thermal=1
			self.thermalLoad_CONV,

			#self.thermalLoad_CHBDY,
			self.chbdye_thermal_load,
			self.chbdyg_thermal_load,
			self.chbdyp_thermal_load,

			#self.thermalLoad_1D,
			self.crod_thermal_load,
			self.cbeam_thermal_load,
			self.ctube_thermal_load,
			self.conrod_thermal_load,
			self.cbar_thermal_load,
			self.cbend_thermal_load,

			#self.thermalLoad_2D_3D,
			self.cquad4_thermal_load,
			self.ctriax6_thermal_load,
			self.cquad8_thermal_load,
			self.ctria3_thermal_load,
			self.ctria6_thermal_load,
			self.ctetra_thermal_load,
			self.cthexa_thermal_load,
			self.cpenta_thermal_load,


			self.thermalLoad_VU,
			self.thermalLoad_VU_3D,
			self.thermalLoad_VUBeam,

			#------------------------------------------

			self.grid_point_stresses, self.grid_point_volume_stresses, self.grid_point_forces,
		]

		for isubcase, res_keys in sorted(iteritems(res_keys_subcase)):
			# print(res_keys)
			for res_key in res_keys:
				if isinstance(res_key, tuple):
					is_compressed = False
				else:
					is_compressed = True

				res_length = self._get_result_length(res_types, res_key)
				if res_length == 0:
					# skipped subcase; no saved results
					continue

				res_format = '*%%-%is SUBCASE=%%i%%s' % res_length
				res_format_vectorized = ' %%-%is SUBCASE=%%i SUBTITLE=%%s %%s' % res_length

				for res_type in res_types:
					if res_key not in res_type:
						continue

					result = res_type[res_key]
					subtitle = result.subtitle
					label = result.label

					header = ['', '']
					header[0] = '	  %-126s\n' % subtitle
					header[1] = '0	 %-32s																	   SUBCASE %-15i\n \n' % (label, isubcase)

					if result.nonlinear_factor is not None:
						header.append('')
					try:
						element_name = ''
						if hasattr(result, 'element_name'):
							element_name = ' - ' + result.element_name

						class_name = result.__class__.__name__
						if hasattr(result, 'data'):
							if not quiet:
								print(res_format_vectorized % (class_name, isubcase, subtitle, element_name))
						else:
							print(res_format % (class_name, isubcase, element_name))

						self.page_num = result.write_f06(f06, header, page_stamp, page_num=self.page_num,
														 is_mag_phase=is_mag_phase, is_sort1=is_sort1)
						assert isinstance(self.page_num, int), 'pageNum=%r' % str(self.page_num)
					except:
						#print("result name = %r" % result.name())
						raise
					if delete_objects:
						del result
					self.page_num += 1

	def oload_resultant_write(self, f, page_stamp, page_num):
		#msg = ''
		#msg += '		*** USER INFORMATION MESSAGE 7310 (VECPRN)\n'
		#msg += '			ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM WILL BE USED AS REFERENCE LOCATION.\n'
		#msg += '			RESULTANTS ABOUT ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM IN SUPERELEMENT BASIC SYSTEM COORDINATES.\n'
		#msg += '	   0												  OLOAD	RESULTANT	   \n'

		"""
		0												  OLOAD	RESULTANT
		  SUBCASE/	LOAD
		  DAREA ID	TYPE	   T1			T2			T3			R1			R2			R3
		0		1	 FX	2.300000E+04	 ----		  ----		  ----	   3.320987E+04 -2.280395E+04
					   FY	   ----	   0.000000E+00	 ----	   0.000000E+00	 ----	   0.000000E+00
					   FZ	   ----		  ----	   0.000000E+00  0.000000E+00  0.000000E+00	 ----
					   MX	   ----		  ----		  ----	   0.000000E+00	 ----		  ----
					   MY	   ----		  ----		  ----		  ----	   0.000000E+00	 ----
					   MZ	   ----		  ----		  ----		  ----		  ----	   0.000000E+00
					 TOTALS  2.300000E+04  0.000000E+00  0.000000E+00  0.000000E+00  3.320987E+04 -2.280395E+04
		#1	MSC.NASTRAN JOB CREATED ON 28-JAN-12 AT 12:52:32					   OCTOBER  22, 2014  MSC.NASTRAN  6/17/05   PAGE	 8
		"""
		isubcase = 1
		Fx = 2.3e4
		msg = '0												  OLOAD	RESULTANT\n'
		msg += '  SUBCASE/	LOAD\n'
		msg += '  DAREA ID	TYPE	   T1			T2			T3			R1			R2			R3\n'
		F = array([[2.3e4, 0., 0.],
				   [0., 0., 0.],
				   [0., 0., 0.]])
		M = array([[0., 3.320987e4, -2.280395e4],
				   [0., 0., 0.],
				   [0., 0., 0.]])
		msg += '0 %8i	 FX	%12.6E	 ----		  ----	   3.320987E+04 -2.280395E+04\n' % (isubcase, Fx)
		msg += '  %8s	 FY	%12s	0.000000E+00	 ----	   0.000000E+00	 ----	   0.000000E+00\n' % ('', '----')
		msg += '  %8s	 FZ	%12s	   ----	   0.000000E+00  0.000000E+00  0.000000E+00	 ----\n' % ('', '----')
		msg += '  %8s	 MX	%12s	   ----		  ----	   0.000000E+00	 ----		  ----\n' % ('', '----')
		msg += '  %8s	 MY	%12s	   ----		  ----		  ----	   0.000000E+00	 ----\n' % ('', '----')
		msg += '  %8s	 MZ	%12s	   ----		  ----		  ----		  ----	   0.000000E+00\n' % ('', '----')
		msg += '  %8s	 TOTALS  2.300000E+04  0.000000E+00  0.000000E+00  0.000000E+00  3.320987E+04 -2.280395E+04\n' % ('')

		msg += page_stamp % page_num
		f.write('\n'.join(msg))
		return page_num + 1


	def _temperature_gradients_and_fluxes(self):
		(subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
		headers = self.skip(2)
		#print "headers = %s" % (headers)
		data = self._read_gradient_fluxes_table()
		return
		if isubcase in self.temperatureGrad:
			self.temperatureGrad[isubcase].addData(data)
		else:
			self.temperatureGrad[isubcase] = TemperatureGradientObject(isubcase, data)
		self.iSubcases.append(isubcase)

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
													 F O R C E S   I N   R O D   E L E M E N T S	 ( C R O D )
		 ELEMENT		   AXIAL									 ELEMENT		   AXIAL
		   ID.			 FORCE		  TORQUE					   ID.			 FORCE		  TORQUE
			 1	   -7.007184E+02   0.0								 2	   -4.900904E+04   0.0
			 3	   -7.141140E+04   0.0
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

								  F O R C E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )		OPTION = BILIN

			ELEMENT					- MEMBRANE  FORCES -					  - BENDING   MOMENTS -			- TRANSVERSE SHEAR FORCES -
			  ID	   GRID-ID	 FX			FY			FXY		   MX			MY			MXY		   QX			QY
				  1	CEN/4  0.0		   0.0		   0.0		  -7.371223E+01 -4.023861E+02 -2.679984E+01  1.315875E+01 -7.356985E+01
						   1  0.0		   0.0		   0.0		  -1.043592E+02 -3.888291E+02 -2.698050E+01  1.315875E+01 -7.356985E+01
						   2  0.0		   0.0		   0.0		  -1.036512E+02 -4.152917E+02 -2.731157E+01  1.315875E+01 -7.356985E+01
						   8  0.0		   0.0		   0.0		  -4.306526E+01 -4.159432E+02 -2.661917E+01  1.315875E+01 -7.356985E+01
						   7  0.0		   0.0		   0.0		  -4.377329E+01 -3.894806E+02 -2.628810E+01  1.315875E+01 -7.356985E+01

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

									   S T R E S S E S   I N   R O D   E L E M E N T S	  ( C R O D )
		  ELEMENT	   AXIAL	   SAFETY	  TORSIONAL	 SAFETY	   ELEMENT	   AXIAL	   SAFETY	  TORSIONAL	 SAFETY
			ID.		STRESS	   MARGIN		STRESS	  MARGIN		 ID.		STRESS	   MARGIN		STRESS	  MARGIN
			   14	2.514247E+04			  1.758725E+02					 15	2.443757E+04			  2.924619E+01
		"""
		(isubcase, transient, dt, data_code) = self._get_rod_header(element_name, element_type, is_strain=False)
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
		(isubcase, transient, dt, data_code) = self._get_rod_header(element_name, element_type, is_strain=True)
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
		* table_code	= 5 (Stress)
		* sort_code	 = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
		* format_code   = 1 (Real)
		* s_code		= 0 (Stress)
		* num_wide	  = 8 (???)
		"""
		(subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
		headers = self.skip(2)
		(stress_bits, s_code) = make_stress_bits(is_strain=is_strain, is_rod_or_solid=True)
		data_code = {
			'analysis_code': analysis_code,
			'device_code': 1, 'table_code': 5, 'sort_code': 0,
			'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
			'stress_bits': stress_bits, 'format_code': 1,
			'element_name': element_name, 'element_type': element_type, 'nonlinear_factor': dt,
			'lsdvmn' : 1,
			'data_names':['lsdvmn']
		}
		return (isubcase, transient, dt, data_code)

	def _read_rod_stress(self):
		"""
		::

									   S T R E S S E S   I N   R O D   E L E M E N T S	  ( C R O D )
		  ELEMENT	   AXIAL	   SAFETY	  TORSIONAL	 SAFETY	   ELEMENT	   AXIAL	   SAFETY	  TORSIONAL	 SAFETY
			ID.		STRESS	   MARGIN		STRESS	  MARGIN		 ID.		STRESS	   MARGIN		STRESS	  MARGIN
			   14	2.514247E+04			  1.758725E+02					 15	2.443757E+04			  2.924619E+01
		"""
		data = []
		while 1:
			line = self.infile.readline()[1:].rstrip('\r\n ')
			sline = [line[0:13], line[13:29], line[29:42], line[42:55], line[55:67], line[67:78], line[78:94], line[94:107], line[107:120], line[120:131]]
			if 'PAGE' in line:
				break
			data_types = [int, float, float, float, float, int, float, float, float, float]
			out = self._parse_line_blanks(sline, data_types)  # line 1
			data.append(out[:5])
			if isinstance(out[5], int):
				data.append(out[5:])
			self.i += 1
		return data

	def _read_spring_stress(self):
		"""
		::

									 S T R A I N S	I N   S C A L A R   S P R I N G S		( C E L A S 2 )
			ELEMENT		 STRAIN		   ELEMENT		 STRAIN		   ELEMENT		 STRAIN		   ELEMENT		 STRAIN
			  ID.							  ID.							  ID.							  ID.
			  20001	  0.0				   20002	  0.0				   20003	  0.0				   20004	  0.0
			  20005	  0.0				   20006	  0.0
		 """
		data = []
		while 1:
			line = self.infile.readline()[1:].rstrip('\r\n ')
			sline = line.strip().split()
			#sline = [line[0:13], line[13:29], line[29:42], line[42:55], line[55:67], line[67:78], line[78:94], line[94:107], line[107:120], line[120:131]]
			if 'PAGE' in line:
				break
			n = len(sline) // 2
			assert len(sline) % 2 == 0, sline

			data_types = [int, float] * n
			out = self._parse_line_blanks(sline, data_types)  # line 1

			while out:
				strain = out.pop()
				eid = out.pop()
				data.append([eid, strain])
			self.i += 1

		return data

	def _stress_in_cbar_elements(self):
		"""
		::

										 S T R E S S E S   I N   B A R   E L E M E N T S		  ( C B A R )
		  ELEMENT		SA1			SA2			SA3			SA4		   AXIAL		  SA-MAX		 SA-MIN	 M.S.-T
			ID.		  SB1			SB2			SB3			SB4		   STRESS		 SB-MAX		 SB-MIN	 M.S.-C
			   12	0.0			0.0			0.0			0.0			1.020730E+04   1.020730E+04   1.020730E+04
					 0.0			0.0			0.0			0.0						   1.020730E+04   1.020730E+04

		* analysis_code = 1 (Statics)
		* device_code   = 1 (Print)
		* table_code	= 5 (Stress)
		* sort_code	 = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
		* format_code   = 1 (Real)
		* s_code		= 0 (Stress)
		* num_wide	  = 8 (???)
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
		(isubcase, transient, dt, data_code) = self._get_bar_header(True)
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

		(stress_bits, s_code) = make_stress_bits(is_strain=is_strain, is_rod_or_solid=True, debug=True)
		data_code = {
			'analysis_code': analysis_code,
			'device_code': 1, 'table_code': 5, 'sort_code': 0,
			'sort_bits': [0, 0, 0], 'num_wide': 8, 's_code': s_code,
			'stress_bits': stress_bits, 'format_code': 1,
			'element_name': 'CBAR', 'element_type': 34,
			'nonlinear_factor': dt,
			'lsdvmn' : 1,
			'data_names':['lsdvmn']
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
			'data_names':['lsdvmn']
		}
		return (isubcase, transient, dt, data_code)

	def _read_bar_stress(self):
		"""
		::

		  ELEMENT		SA1			SA2			SA3			SA4		   AXIAL		  SA-MAX		 SA-MIN	 M.S.-T
			ID.		  SB1			SB2			SB3			SB4		   STRESS		 SB-MAX		 SB-MIN	 M.S.-C
			   12	0.0			0.0			0.0			0.0			1.020730E+04   1.020730E+04   1.020730E+04
					 0.0			0.0			0.0			0.0						   1.020730E+04   1.020730E+04
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
		  ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS	INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)	  MAX
			ID	  ID	NORMAL-1	 NORMAL-2	 SHEAR-12	 SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE	MAJOR		MINOR		SHEAR
			  181	1   3.18013E+04  5.33449E+05  1.01480E+03   -7.06668E+01  1.90232E+04   89.88  5.33451E+05  3.17993E+04  2.50826E+05
			  181	2   1.41820E+05  1.40805E+05  1.25412E+05   -1.06000E+02  2.85348E+04   44.88  2.66726E+05  1.58996E+04  1.25413E+05

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
		  ELEMENT  PLY  STRESSES IN FIBER AND MATRIX DIRECTIONS	INTER-LAMINAR  STRESSES  PRINCIPAL STRESSES (ZERO SHEAR)	  MAX
			ID	  ID	NORMAL-1	 NORMAL-2	 SHEAR-12	 SHEAR XZ-MAT  SHEAR YZ-MAT  ANGLE	MAJOR		MINOR		SHEAR
			  181	1   3.18013E+04  5.33449E+05  1.01480E+03   -7.06668E+01  1.90232E+04   89.88  5.33451E+05  3.17993E+04  2.50826E+05
			  181	2   1.41820E+05  1.40805E+05  1.25412E+05   -1.06000E+02  2.85348E+04   44.88  2.66726E+05  1.58996E+04  1.25413E+05

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
			'data_names':['lsdvmn'],
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
		  ELEMENT	  FIBER			   STRESSES IN ELEMENT COORD SYSTEM			 PRINCIPAL STRESSES (ZERO SHEAR)
			ID.	   DISTANCE		   NORMAL-X	   NORMAL-Y	  SHEAR-XY	   ANGLE		 MAJOR		   MINOR		VON MISES
				8   -1.250000E-01	 -1.303003E+02   1.042750E+04  -1.456123E+02   -89.2100	1.042951E+04   -1.323082E+02   1.049629E+04
					 1.250000E-01	 -5.049646E+02   1.005266E+04  -2.132942E+02   -88.8431	1.005697E+04   -5.092719E+02   1.032103E+04

		* analysis_code = 1 (Statics)
		* device_code   = 1 (Print)
		* table_code	= 5 (Stress)
		* sort_code	 = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
		* format_code   = 1 (Real)
		* s_code		= 0 (Stress)
		* num_wide	  = 8 (???)
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
		* table_code	= 5 (Stress)
		* sort_code	 = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
		* format_code   = 1 (Real)
		* s_code		= 0 (Stress)
		* num_wide	  = 8 (???)
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
			'data_names':['lsdvmn'],
			'lsdvmn': 1,
			}
		if transient is not None:
			data_code['name'] = transient[0]
		return (isubcase, transient, data_code)

	def _read_tri_stress(self, eType):
		"""
		::

				  ELEMENT	  FIBER			   STRESSES IN ELEMENT COORD SYSTEM			 PRINCIPAL STRESSES (ZERO SHEAR)
					ID.	   DISTANCE		   NORMAL-X	   NORMAL-Y	  SHEAR-XY	   ANGLE		 MAJOR		   MINOR		VON MISES
						8   -1.250000E-01	 -1.303003E+02   1.042750E+04  -1.456123E+02   -89.2100	1.042951E+04   -1.323082E+02   1.049629E+04
							 1.250000E-01	 -5.049646E+02   1.005266E+04  -2.132942E+02   -88.8431	1.005697E+04   -5.092719E+02   1.032103E+04
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

							   S T R E S S E S   I N   Q U A D R I L A T E R A L   E L E M E N T S   ( Q U A D 4 )		OPTION = BILIN

		  ELEMENT			  FIBER			STRESSES IN ELEMENT COORD SYSTEM		 PRINCIPAL STRESSES (ZERO SHEAR)
			ID	  GRID-ID   DISTANCE		NORMAL-X	  NORMAL-Y	  SHEAR-XY	  ANGLE		MAJOR		 MINOR	   VON MISES
				6	CEN/4  -1.250000E-01  -4.278394E+02  8.021165E+03 -1.550089E+02   -88.9493   8.024007E+03 -4.306823E+02  8.247786E+03
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
			is_von_mises = not is_max_shear
			assert result.is_max_shear() == is_max_shear
			assert result.is_von_mises() == is_von_mises
		self.iSubcases.append(isubcase)

	def _get_quad_header(self, nheader_lines, element_name, element_type, is_strain):
		(subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
		headers = self.skip(nheader_lines)
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
			'data_names':['lsdvmn'],
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

		data = self._read_3d_stress(element_name, nnodes, is_strain=False)
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

		data = self._read_3d_stress(element_name, nnodes, is_strain=True)
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
		* table_code	= 5 (Stress/Strain)
		* sort_code	 = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
		* format_code   = 1 (Real)
		* s_code		= 0 (Stress/Strain)
		* num_wide	  = 8 (???)
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
			'data_names':['lsdvmn']
		}
		if transient is not None:
			data_code['name'] = transient[0]
		return (isubcase, transient, data_code)

	def _read_3d_stress(self, etype, nnodes, is_strain=False):
		data = []
		etype2 = etype + str(nnodes)
		while 1:
			line = self.infile.readline().rstrip('\n\r')  # [1:]
					#			  CENTER		 X		  #		  XY			 #		A		 #
			sline = [line[1:17], line[17:23], line[23:28], line[28:43], line[43:47], line[47:63], line[63:66], line[66:80], line[80:83], line[83:88], line[88:93], line[93:98], line[99:113], line[113:130]]
			sline = [s.strip() for s in sline]
			if 'PAGE' in line:
				break
			elif '' is not sline[0]:
				sline = [etype2] + sline
			data.append(sline)

		return data

def make_stress_bits(is_fiber_distance=False, is_max_shear=True, is_strain=True, is_rod_or_solid=False, debug=False):
	"""
	Therefore, stress_code can be one of the following values:
	+------+---------+----------------------------------------------+
	|Value | On bits | Description								  |
	+------+---------+----------------------------------------------+
	|  0   | 0 0 0 0 | Stress maximum shear or octahedral		   |
	|  1   | 0 0 0 1 | Stress von Mises							 |
	|  10  | 1 0 1 0 | Strain Curvature maximum shear or octahedral |
	|  11  | 1 0 1 1 | Strain Curvature von Mises				   |
	|  14  | 1 1 1 0 | Strain Fibre maimum shear or octahedral	  |
	|  15  | 1 1 1 1 | Strain Fibre von Mises					   |
	+------+---------+----------------------------------------------+
	"""
	#if debug:
	msg = ("  is_max_shear=%s is_fiber_distance=%s is_strain=%s is_rod_or_solid=%s"
		   % (is_max_shear, is_fiber_distance, is_strain, is_rod_or_solid))
	#print(msg)
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
		#(is_max_shear, is_fiber_distance, is_strain, is_rod_or_solid), scode
		(True,  False, True,  True) : ([0, 1, 1, 1, 0], None), # bar (rod/solid) strain
		(False, False, True,  True) : ([0, 1, 1, 1, 1], None), # solid (rod/solid) strain
		(True,  False, False, True) : ([0, 0, 1, 0, 0], None), # rod (crod/solid) stress
		(False, False, False, True) : ([0, 0, 1, 0, 1], None), # solid (crod/solid) stress

		#(True,  False, False, True) : ([0, 0, 0, 0, 0], 0),  # 0,  rod/csolid
		#(False, False, False, True) : ([0, 0, 0, 0, 1], 1),  # 1,  rod/csolid strain
		#(False, False, True,  True) : ([0, 1, 0, 1, 1], 1),   # ???, rod/csolid strain
		#(True,  False, True,  True) : ([0, 1, 0, 1, 1], 1), # cbar strain...probably not correct...

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
		raise RuntimeError(msg)

	if is_strain:
		# strain
		assert stress_bits[1] == 1, stress_bits
	else:
		# stress
		assert stress_bits[1] == 0, stress_bits
	assert stress_bits[1] == stress_bits[3], stress_bits

	# backwards
	if is_fiber_distance:
		# fiber distance
		assert stress_bits[2] == 0, stress_bits
	else:
		# curvature
		assert stress_bits[2] == 1, stress_bits

	# backwards
	if is_max_shear:
		assert stress_bits[4] == 0, stress_bits
	else:
		assert stress_bits[4] == 1, stress_bits

	# def is_curvature(self):
		# return True if self.stress_bits[2] == 1 else False
	# def is_max_shear(self):
		# return True if self.stress_bits[4] == 0 else False

	# def is_fiber_distance(self):
		# return not self.is_curvature()
		# return True if self.stress_bits[2] == 0 else False

	# def is_von_mises(self):
		# return not self.is_max_shear()
		# return True if self.stress_bits[4] == 1 else False


	#if is_max_shear==False:
	#	stress_bits[4] = 1 # Von Mises
	#if is_strain:
	#	#stress_bits[1] = stress_bits[3] = 1 # Strain
	#	stress_bits[1] = stress_bits[3] = 1 # Strain
	#if is_fiber_distance:
	#	stress_bits[2] = 1 # FiberDistance
	#print stress_bits
	#s_code = 0
	#for i,bit in enumerate(stress_bits):
	#	s_code += bit*2**i

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

	def _read_f06_table(self, data_types, debug=False):
		pass
	def __init__(self):
		self.iSubcases = []
		self.i = 0

	def _load_vector(self):
		(subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
		headers = self.skip(2)
		#print "headers = %s" % (headers)

		data = self._real_f06_real_table_data(allow_blanks=False)

		data_code = {
			'analysis_code': analysis_code,
			'device_code': 1,

			'table_code': 1, # ???
			'table_name': 'OPG',
			'format_code': 1, # ???

			'sort_code': 0,
			'sort_bits': [0, 0, 0], 'num_wide': 8,
			'nonlinear_factor': dt,
			'data_names':['lsdvmn'],
			'lsdvmn': 1,
			}

		if isubcase in self.load_vectors:
			self.load_vectors[isubcase].add_f06_data(data, transient)
		else:
			is_sort1 = True
			spc = RealLoadVector(data_code, is_sort1, isubcase, dt)
			spc.add_f06_data(data, transient)
			self.load_vectors[isubcase] = spc
		self.iSubcases.append(isubcase)

	def _forces_of_single_point_constraints(self):
		(subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
		headers = self.skip(2)
		#print "headers = %s" %(headers)

		data = self._real_f06_real_table_data(allow_blanks=False)

		data_code = {
			'analysis_code': analysis_code,
			'device_code': 1, 'table_code': 3, 'sort_code': 0,
			'sort_bits': [0, 0, 0], 'num_wide': 8, 'table_name': 'OQG',
			'nonlinear_factor': dt,
			'format_code': 3,  # ???
			'data_names':['lsdvmn'],
			'lsdvmn': 1,
		}

		if isubcase in self.spc_forces:
			self.spc_forces[isubcase].add_f06_data(data, transient)
		else:
			is_sort1 = True
			spc = RealSPCForces(data_code, is_sort1, isubcase, dt)
			spc.add_f06_data(data, transient)
			self.spc_forces[isubcase] = spc
		self.iSubcases.append(isubcase)

	def _forces_of_multi_point_constraints(self):
		(subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
		headers = self.skip(2)
		#print "headers = %s" %(headers)

		data = self._real_f06_real_table_data(allow_blanks=False)

		data_code = {
			'analysis_code': analysis_code,
			'device_code': 1, 'table_code': 39,
			'sort_code': 0, 'sort_bits': [0, 0, 0], 'num_wide': 8,
			'table_name': 'OQG', 'nonlinear_factor': dt,
			'data_names':['lsdvmn']
		}

		if isubcase in self.mpc_forces:
			self.mpc_forces[isubcase].add_f06_data(data, transient)
		else:
			is_sort1 = True
			mpc = RealMPCForces(data_code, is_sort1, isubcase, dt)
			mpc.add_f06_data(data, transient)
			self.mpc_forces[isubcase] = mpc
		self.iSubcases.append(isubcase)


	def _read_f06_table(self, data_types, debug=False):
		pass

	def _real_eigenvectors(self, marker):
		"""
		Reads real eigenvector table accounting for blank entries

		::
																												 SUBCASE 1
		  EIGENVALUE =  6.158494E+07
			  CYCLES =  1.248985E+03		 R E A L   E I G E N V E C T O R   N O .		  1

		  POINT ID.   TYPE		  T1			 T2			 T3			 R1			 R2			 R3
				 1	  G	  2.547245E-17  -6.388945E-16   2.292728E+00  -1.076928E-15   2.579163E-17   0.0
			  2002	  G	 -6.382321E-17  -1.556607E-15   3.242408E+00  -6.530917E-16   1.747180E-17   0.0
			  2003	  G	 -6.382321E-17  -1.556607E-15   3.242408E+00
			  2004	  S	 -6.382321E-17  -1.556607E-15   3.242408E+00

		* analysis_code = 2 (Normal modes)
		* table_code	= 7 (Eigenvector)
		* device_code   = 1 (Print)
		* sort_code	 = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
		* format_code   = 1 (Real)
		* #s_code		= 0 (Stress)
		* num_wide	  = 8 (???)
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
												 C O M P L E X   E I G E N V E C T O R   NO.		  5
																	(REAL/IMAGINARY)

				POINT ID.   TYPE		  T1			 T2			 T3			 R1			 R2			 R3
		  0			1	  G	  0.0			0.0			0.0		   -3.041936E-02  -1.922321E-01   0.0
									 0.0			0.0			0.0			0.0			0.0			0.0
		  0		  227	  S	  1.418276E-03  -2.095675E-04  -1.663478E-03   2.633889E-03   4.171373E-14  -2.633889E-03
									 0.0			0.0			0.0			0.0			0.0			0.0
		  0		  233	  S	  1.663478E-03   2.095675E-04  -1.418276E-03
									 0.0			0.0			0.0

		* analysis_code = 2 (Normal modes)
		* table_code	= 7 (Eigenvector)
		* device_code   = 1 (Print)
		* sort_code	 = ? (Sort2,Real/Imag or Mag/Phase,Sorted Results) => sort_bits = [0,0,0]
		* format_code   = 1 (???)
		* num_wide	  = 16 (8*2)
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
		num_fields = 6	 # the number of fields (T1, T2, T3, R1, R2, R3)
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
		num_fields = 6	 # the number of fields (T1, T2, T3, R1, R2, R3)
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

		  POINT ID.   TYPE		  T1			 T2			 T3			 R1			 R2			 R3
				 1	  G	  9.663032E-05   0.0		   -2.199001E-04   0.0		   -9.121119E-05   0.0
				 2	  G	  0.0			0.0			0.0			0.0			0.0			0.0
				 3	  G	  0.0			0.0			0.0			0.0			0.0			0.0

		* analysis_code = 1 (Statics)
		* device_code   = 1 (Print)
		* table_code	= 1 (Displacement)
		* sort_code	 = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
		* num_wide	  = 8 (???)
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
			POINT-ID =	   101
											 C O M P L E X   D I S P L A C E M E N T   V E C T O R
															   (MAGNITUDE/PHASE)

			FREQUENCY   TYPE		  T1			 T2			 T3			 R1			 R2			 R3
		  2.000000E+01	 G	  3.242295E-16   1.630439E-01   1.630439E-01   1.691497E-17   1.362718E-01   1.362718E-01
								  196.0668		90.0000	   180.0000		63.4349	   180.0000	   270.0000

		* table_code	= 1 (Displacement)
		* format_code   = 3 (Magnitude/Phase)
		* sort_bits	 = [0,1,1]  (Sort1,Real/Imaginary,RandomResponse)
		* analysis_code = 5 (Frequency)
		* sort_code	 = 2 (Random Response)
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

		  POINT ID.   TYPE	  ID   VALUE	 ID+1 VALUE	 ID+2 VALUE	 ID+3 VALUE	 ID+4 VALUE	 ID+5 VALUE
				 1	  S	  1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03
				 7	  S	  1.300000E+03   1.300000E+03   1.300000E+03   1.300000E+03
		  analysis_code = 1 (Statics)
		  device_code   = 1 (Print)
		  table_code	= 1 (Displacement/Temperature)
		  sort_code	 = 0 (Sort2,Real,Sorted Results) => sort_bits = [0,0,0]
		  format_code   = 1 (Real)
		  s_code		= 0 (Stress)
		  num_wide	  = 8 (???)
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


	def _real_eigenvalues(self):
		"""
		::

													 R E A L   E I G E N V A L U E S
		   MODE	EXTRACTION	  EIGENVALUE			RADIANS			 CYCLES			GENERALIZED		 GENERALIZED
			NO.	   ORDER																	   MASS			  STIFFNESS
				1		 1		6.158494E+07		7.847607E+03		1.248985E+03		1.000000E+00		6.158494E+07
		"""
		self.title = None
		(subcase_name, isubcase, transient, dt, analysis_code, is_sort1) = self._read_f06_subcase_header()
		Title = None
		line1 = self.infile.readline().strip(); self.i += 1
		if line1 != 'MODE	EXTRACTION	  EIGENVALUE			RADIANS			 CYCLES			GENERALIZED		 GENERALIZED':
			Title = line1
			line1 = self.infile.readline().strip(); self.i += 1
		line2 = self.infile.readline().strip(); self.i += 1

		#MODE	EXTRACTION	  EIGENVALUE			RADIANS			 CYCLES			GENERALIZED		 GENERALIZED
		# NO.	   ORDER																	   MASS			  STIFFNESS
		#	 1		 1		1.018377E-03		3.191203E-02		5.078956E-03		1.000000E+00		1.018377E-03
		#print(line1)
		#print(line2)
		#headers = self.skip(2)
		#print(headers)
		data = self._read_f06_table([int, int, float, float, float, float, float])

		self.eigenvalues[self.title] = RealEigenvalues(Title)
		self.eigenvalues[self.title].add_f06_data(data)

	def _complex_eigenvalue_summary(self):
		"""
		::

								 C O M P L E X   E I G E N V A L U E   S U M M A R Y
		  ROOT	 EXTRACTION				  EIGENVALUE					 FREQUENCY			  DAMPING
		   NO.		ORDER			 (REAL)		   (IMAG)				(CYCLES)			COEFFICIENT
			   1		   6		  0.0			  6.324555E+01		  1.006584E+01		  0.0
			   2		   5		  0.0			  6.324555E+01		  1.006584E+01		  0.0
		"""
		#(subcaseName,isubcase,transient,dt,analysis_code,is_sort1) = self.readSubcaseNameID()
		isubcase = 1  # .. todo:: fix this...

		headers = self.skip(2)
		data = self._read_f06_table([int, int, float, float, float, float])

		if self.title in self.eigenvalues:
			self.eigenvalues[self.title].add_f06_data(data)
		else:
			self.eigenvalues[self.title] = ComplexEigenvalues(self.title)
			self.eigenvalues[self.title].add_f06_data(data)
		self.iSubcases.append(isubcase)
