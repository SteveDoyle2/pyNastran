from __future__ import print_function
from six import iteritems, string_types, binary_type, text_type
from collections import defaultdict
from numpy import unique, int32, int64

from pyNastran import is_release
from pyNastran.op2.tables.grid_point_weight import GridPointWeight
from pyNastran.f06.f06_formatting import get_key0
from pyNastran.utils import object_attributes, integer_types
from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.bdf.case_control_deck import CaseControlDeck

try:
    import pandas as pd
except ImportError:
    pass


class OP2_F06_Common(object):
    def __init__(self):
        #: a dictionary that maps an integer of the subcaseName to the
        #: subcaseID
        self.iSubcaseNameMap = {}
        self.subtitles = defaultdict(list)
        self.case_control_deck = CaseControlDeck([], log=self.log)
        self.labels = {}

        self.make_geom = False

        #: BDF Title
        self.title = None

        self.page_num = 1

        self.iSubcases = []
        self.__objects_vector_init__()
        self.__objects_init__()
        self.__objects_common_init__()

    def deprecated(self, old_name, new_name, deprecated_version):
        """allows for simple OP2 vectorization"""
        return deprecated(old_name, new_name, deprecated_version, levels=[0, 1, 2])

    def __objects_vector_init__(self):
        """
        All OUG table is simple to vectorize, so we declere it in __objects_init__
        On the other hand, the rodForces object contains CROD/CTUBE/CONROD elements.
        It is difficult to handle initializing the CRODs/CONRODs given a
        mixed type case, so we split out the elements.
        """
        #======================================================================
        # rods
        self.crod_force = {}
        self.conrod_force = {}
        self.ctube_force = {}

        self.crod_stress = {}
        self.conrod_stress = {}
        self.ctube_stress = {}

        self.crod_strain = {}
        self.conrod_strain = {}
        self.ctube_strain = {}

        #======================================================================
        # springs
        self.celas1_force = {}
        self.celas2_force = {}
        self.celas3_force = {}
        self.celas4_force = {}

        self.celas1_stress = {}
        self.celas2_stress = {}
        self.celas3_stress = {}
        self.celas4_stress = {}

        self.celas1_strain = {}
        self.celas2_strain = {}
        self.celas3_strain = {}
        self.celas4_strain = {}

        #======================================================================
        self.ctetra_stress = {}
        self.cpenta_stress = {}
        self.chexa_stress = {}

        self.ctetra_strain = {}
        self.cpenta_strain = {}
        self.chexa_strain = {}
        #======================================================================

        # bars/beams
        self.cbar_force = {}
        self.cbar_stress = {}
        self.cbar_strain = {}

        self.cbar_force_10nodes = {}
        self.cbar_stress_10nodes = {}
        self.cbar_strain_10nodes = {}

        self.cbeam_force = {}
        self.cbeam_stress = {}
        self.cbeam_strain = {}

        #======================================================================
        # shells
        self.ctria3_force = {}
        self.ctria6_force = {}
        self.ctriar_force = {}

        self.cquad4_force = {}
        self.cquad8_force = {}
        self.cquadr_force = {}

        self.ctria3_stress = {}
        self.ctria6_stress = {}
        self.cquad4_stress = {}
        self.cquad8_stress = {}
        self.cquadr_stress = {}
        self.ctriar_stress = {}

        self.ctria3_strain = {}
        self.ctria6_strain = {}
        self.cquad4_strain = {}
        self.cquad8_strain = {}
        self.cquadr_strain = {}
        self.ctriar_strain = {}

        self.cquad4_composite_stress = {}
        self.cquad8_composite_stress = {}
        self.cquadr_composite_stress = {}
        self.ctria3_composite_stress = {}
        self.ctria6_composite_stress = {}
        self.ctriar_composite_stress = {}

        self.cquad4_composite_strain = {}
        self.cquad8_composite_strain = {}
        self.cquadr_composite_strain = {}
        self.ctria3_composite_strain = {}
        self.ctria6_composite_strain = {}
        self.ctriar_composite_strain = {}

        self.cplstn3_stress = {}
        self.cplstn4_stress = {}
        self.cplstn6_stress = {}
        self.cplstn8_stress = {}
        self.cplsts3_stress = {}
        self.cplsts4_stress = {}
        self.cplsts6_stress = {}
        self.cplsts8_stress = {}

        self.cplstn3_strain = {}
        self.cplstn4_strain = {}
        self.cplstn6_strain = {}
        self.cplstn8_strain = {}
        self.cplsts3_strain = {}
        self.cplsts4_strain = {}
        self.cplsts6_strain = {}
        self.cplsts8_strain = {}

        self.cshear_stress = {}
        self.cshear_strain = {}
        self.cshear_force = {}

        #: OES - CBEAM 94
        self.nonlinear_cbeam_stress = {}

        # bushing
        self.cbush_stress = {}
        self.cbush_strain = {}
        self.nonlinear_cbush_stress = {}  # CBUSH 226
        self.cbush1d_stress_strain = {}

        #======================================================================

    def __objects_common_init__(self):
        #: the date the job was run on
        self.date = None

        # SOL 200
        self.convergence_data = None
        self.weight_response = None
        self.stress_response = None
        self.strain_response = None
        self.composite_stress_response = None
        self.composite_strain_response = None
        self.flutter_response = None

        #: Grid Point Weight Table
        #: create with:
        #:   PARAM   GRDPNT    0  (required for F06/OP2)
        #:   PARAM   POSTEXT YES  (required for OP2)
        self.grid_point_weight = GridPointWeight()
        self.oload_resultant = None

        #: ESE
        self.eigenvalues = {}

        #self.convergence_history = {}
        #self.response1_table = {}

    def __objects_init__(self):
        """More variable declarations"""
        #: the date the job was run on
        self.date = None

        #: Grid Point Weight Table
        #: create with:
        #:   PARAM   GRDPNT    0  (required for F06/OP2)
        #:   PARAM   POSTEXT YES  (required for OP2)
        self.grid_point_weight = GridPointWeight()

        #: ESE
        self.eigenvalues = {}

        #: OUG - displacement
        self.displacements = {}           # tCode=1 thermal=0
        self.displacementsPSD = {}        # random
        self.displacementsATO = {}        # random
        self.displacementsRMS = {}        # random
        self.displacementsCRM = {}        # random
        self.displacementsNO = {}         # random
        self.displacements_scaled = {}    # tCode=1 thermal=8

        #: OUP

        self.displacement_scaled_response_spectra_NRL = {}  # thermal=8
        self.displacement_scaled_response_spectra_ABS = {}  # thermal=2
        self.displacement_scaled_response_spectra_SRSS = {} # thermal=4
        #self.displacement_scaled_response_spectra_PSD = {}
        #self.displacement_scaled_response_spectra_ATO = {}
        #self.displacement_scaled_response_spectra_RMS = {}
        #self.displacement_scaled_response_spectra_CRM = {}
        #self.displacement_scaled_response_spectra_NO = {}


        #: OUG - velocity
        self.velocities = {}              # tCode=10 thermal=0
        self.velocitiesPSD = {}
        #self.velocity_scaled_response_spectra_NRL = {}
        self.velocity_scaled_response_spectra_ABS = {}
        #self.velocity_scaled_response_spectra_PSD = {}
        #self.velocity_scaled_response_spectra_ATO = {}
        #self.velocity_scaled_response_spectra_RMS = {}
        #self.velocity_scaled_response_spectra_CRM = {}
        #self.velocity_scaled_response_spectra_NO = {}

        #: OUG - acceleration
        self.accelerations = {}            # tCode=11 thermal=0
        self.accelerationsPSD = {}
        self.acceleration_scaled_response_spectra_NRL = {}
        self.acceleration_scaled_response_spectra_ABS = {}
        #self.acceleration_scaled_response_spectra_PSD = {}
        #self.acceleration_scaled_response_spectra_ATO = {}
        #self.acceleration_scaled_response_spectra_RMS = {}
        #self.acceleration_scaled_response_spectra_CRM = {}
        #self.acceleration_scaled_response_spectra_NO = {}


        #: OUG - temperatures
        self.temperatures = {}           # tCode=1 thermal=1

        #: OUG - eigenvectors
        self.eigenvectors = {}            # tCode=7 thermal=0
        self.eigenvectors_RADCONS = {}
        self.eigenvectors_RADEFFM = {}
        self.eigenvectors_RADEATC = {}

        # OEF - Forces - tCode=4 thermal=0

        self.cbend_force = {}
        self.cbush_force = {}
        self.coneax_force = {}

        self.cdamp1_force = {}
        self.cdamp2_force = {}
        self.cdamp3_force = {}
        self.cdamp4_force = {}

        self.celas1_force = {}
        self.celas2_force = {}
        self.celas3_force = {}
        self.celas4_force = {}

        self.cgap_force = {}

        #self.solidPressureForces = {}
        self.chexa_pressure_force = {}
        self.cpenta_pressure_force = {}
        self.ctetra_pressure_force = {}

        self.cvisc_force = {}

        self.force_VU = {}
        self.force_VU_2D = {}

        #OEF - Fluxes - tCode=4 thermal=1
        self.thermalLoad_CONV = {}

        #self.thermalLoad_CHBDY = {}
        self.chbdye_thermal_load = {}
        self.chbdyg_thermal_load = {}
        self.chbdyp_thermal_load = {}

        #self.thermalLoad_1D = {}
        self.crod_thermal_load = {}
        self.cbeam_thermal_load = {}
        self.ctube_thermal_load = {}
        self.conrod_thermal_load = {}
        self.cbar_thermal_load = {}
        self.cbend_thermal_load = {}

        #self.thermalLoad_2D_3D = {}
        self.cquad4_thermal_load = {}
        self.ctriax6_thermal_load = {}
        self.cquad8_thermal_load = {}
        self.ctria3_thermal_load = {}
        self.ctria6_thermal_load = {}
        self.ctetra_thermal_load = {}
        self.cthexa_thermal_load = {}
        self.cpenta_thermal_load = {}

        self.thermalLoad_VU = {}
        self.thermalLoad_VU_3D = {}
        self.thermalLoad_VUBeam = {}
        #self.temperatureForces = {}

        # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)

        #: OES - CTRIAX6
        self.ctriax_stress = {}
        self.ctriax_strain = {}

        #: OES - nonlinear CROD/CONROD/CTUBE stress/strain
        self.nonlinear_crod_stress = {}
        self.nonlinear_crod_strain = {}

        self.nonlinear_ctube_stress = {}
        self.nonlinear_ctube_strain = {}

        self.nonlinear_conrod_stress = {}
        self.nonlinear_conrod_strain = {}

        #: OESNLXR - CTRIA3/CQUAD4 stress
        #self.nonlinearPlateStress = {}
        #: OESNLXR - CTRIA3/CQUAD4 strain
        #self.nonlinearPlateStrain = {}
        #self.hyperelastic_plate_stress = {}
        self.hyperelastic_cquad4_strain = {}

        self.nonlinear_cquad4_stress = {}
        self.nonlinear_ctria3_stress = {}

        self.nonlinear_cquad4_strain = {}
        self.nonlinear_ctria3_strain = {}


        #: OES - CELAS1 224, CELAS3 225,
        self.nonlinear_celas1_stress = {}
        self.nonlinear_celas3_stress = {}

        #: OES - GAPNL 86
        self.nonlinear_cgap_stress = {}

        # OQG - spc/mpc forces
        self.spc_forces = {}  # tCode=3?
        self.spc_forces_scaled_response_spectra_NRL = {}
        self.spc_forcesPSD = {}
        self.spc_forcesATO = {}
        self.spc_forcesRMS = {}

        self.mpc_forces = {}  # tCode=39
        self.mpc_forcesPSD = {}
        self.mpc_forcesATO = {}
        self.mpc_forcesRMS = {}
        self.mpc_forces_RAQCONS = {}
        #self.mpc_forces_RAQEATC = {}

        # OQG - thermal forces
        self.thermal_gradient_and_flux = {}

        #: OGF - grid point forces
        self.grid_point_forces = {}  # tCode=19

        #: OGS1 - grid point stresses
        self.grid_point_stresses = {}       # tCode=26
        self.grid_point_volume_stresses = {}  # tCode=27

        #: OPG - summation of loads for each element
        self.load_vectors = {}       # tCode=2  thermal=0
        self.thermal_load_vectors = {}  # tCode=2  thermal=1
        self.applied_loads = {}       # tCode=19 thermal=0
        self.force_vectors = {}       # tCode=12 thermal=0

        #: OEE - strain energy density
        #self.strain_energy = {}  # tCode=18
        self.cquad4_strain_energy = {}
        self.cquad8_strain_energy = {}
        self.cquadr_strain_energy = {}
        self.cquadx_strain_energy = {}

        self.ctria3_strain_energy = {}
        self.ctria6_strain_energy = {}
        self.ctriar_strain_energy = {}
        self.ctriax_strain_energy = {}
        self.ctriax6_strain_energy = {}

        self.ctetra_strain_energy = {}
        self.cpenta_strain_energy = {}
        self.chexa_strain_energy = {}
        self.cpyram_strain_energy = {}

        self.crod_strain_energy = {}
        self.ctube_strain_energy = {}
        self.conrod_strain_energy = {}

        self.cbar_strain_energy = {}
        self.cbeam_strain_energy = {}

        self.cgap_strain_energy = {}
        self.celas1_strain_energy = {}
        self.celas2_strain_energy = {}
        self.celas3_strain_energy = {}
        self.celas4_strain_energy = {}
        self.cdum8_strain_energy = {}
        self.cbush_strain_energy = {}
        #self.chexa8fd_strain_energy = {}
        self.cbend_strain_energy = {}
        self.dmig_strain_energy = {}
        self.genel_strain_energy = {}
        self.cshear_strain_energy = {}

    def _get_result_length(self, res_types, res_key):
        """
        gets the length of the output data so we can line up:

          RealCRodStrain  - CROD
          RealCTubeStrain - CTUBE
        """
        res_length = 0
        for res_type in res_types:
            if not res_type:
                continue
            key0 = next(iter(res_type))
            if not isinstance(key0, integer_types) and not isinstance(res_key, integer_types):
                if not type(key0) == type(res_key):
                    msg = 'bad compression check...keys0=%s type(key0)=%s res_key=%s type(res_key)=%s' % (
                        key0, type(key0), res_key, type(res_key))
                    raise RuntimeError(msg)

            #print('res_type.keys()=%s' % res_type.keys())
            # res_key_list = res_key[:-1] + [res_key[-1]]
            # res_key = tuple(res_key_list)

            if res_key in res_type:
                # the res_key is
                result = res_type[res_key]
                class_name = result.__class__.__name__
                res_length = max(len(class_name), res_length)
                #print('continue')
                #break
                continue
            elif len(res_type) != 0:
                #print('  not valid')
                # get the 0th key in the dictionary, where key0 is arbitrary
                key0 = get_key0(res_type)
                #print('  key0 = ', key0)

                # extract displacement[0]
                result = res_type[key0]

                # get the class name
                class_name = result.__class__.__name__
                res_length = max(len(class_name), res_length)

                if not is_release:
                    print('%s - results not found...key=%s' % (class_name, res_key))
            else:  # empty result
                #print('else')
                pass
        #print('res_length =', res_length)
        return res_length

    def get_table_types(self):
        """
        Gets the names of the results.
        """
        table_types = [
            # OUG - displacement
            'displacements',
            'displacementsPSD',
            'displacementsATO',
            'displacementsRMS',
            'displacementsCRM',
            'displacementsNO',
            'displacements_scaled',

            # OUG - temperatures
            'temperatures',

            # OUG - eigenvectors
            'eigenvectors',
            'eigenvectors_RADCONS',
            'eigenvectors_RADEFFM',
            'eigenvectors_RADEATC',


            # OUG - velocity
            'velocities',
            'velocitiesPSD',

            # OUG - acceleration
            'accelerations',
            'accelerationsPSD',

            # OQG - spc/mpc forces
            'spc_forces', 'spc_forcesPSD', 'spc_forcesATO', 'spc_forcesRMS',
            'spc_forces_scaled_response_spectra_NRL',

            'mpc_forces', 'mpc_forcesPSD', 'mpc_forcesATO', 'mpc_forcesRMS',
            'mpc_forces_RAQCONS',
            #'mpc_forces_RAQEATC',

            'thermal_gradient_and_flux',

            # OGF - grid point forces
            'grid_point_forces',

            # OPG - summation of loads for each element
            'load_vectors',
            'thermal_load_vectors',
            'applied_loads',
            'force_vectors',

            # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
            # OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
            'celas1_stress',
            'celas2_stress',
            'celas3_stress',
            'celas4_stress',

            # OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
            'celas1_strain',
            'celas2_strain',
            'celas3_strain',
            'celas4_strain',

            # OES - isotropic CROD/CONROD/CTUBE stress
            'crod_stress',
            'conrod_stress',
            'ctube_stress',

            # OES - isotropic CROD/CONROD/CTUBE strain
            'crod_strain',
            'conrod_strain',
            'ctube_strain',

            # OES - isotropic CBAR stress/strain
            'cbar_stress',
            'cbar_strain',
            'cbar_force',

            'cbar_stress_10nodes',
            'cbar_strain_10nodes',
            'cbar_force_10nodes',

            # OES - isotropic CBEAM stress/strain
            'cbeam_stress',
            'cbeam_strain',
            'cbeam_force',
            'nonlinear_cbeam_stress',
            #'nonlinear_cbeam_strain',


            # OES - isotropic CTRIA3/CQUAD4 stress
            'ctria3_stress',
            'ctriar_stress',
            'ctria6_stress',

            'cquadr_stress',
            'cquad4_stress',
            'cquad8_stress',

            # OES - isotropic CTRIA3/CQUAD4 strain
            'ctria3_strain',
            'ctriar_strain',
            'ctria6_strain',

            'cquadr_strain',
            'cquad4_strain',
            'cquad8_strain',


            # OES - isotropic CTETRA/CHEXA/CPENTA stress
            'ctetra_stress',
            'chexa_stress',
            'cpenta_stress',

            # OES - isotropic CTETRA/CHEXA/CPENTA strain
            'ctetra_strain',
            'chexa_strain',
            'cpenta_strain',

            # OES - CSHEAR stress/strain
            'cshear_stress',
            'cshear_strain',

            # OES - GAPNL 86
            'nonlinear_cgap_stress',
            # OES - CBUSH 226
            'nonlinear_cbush_stress',

            'cplstn3_stress',
            'cplstn4_stress',
            'cplstn6_stress',
            'cplstn8_stress',
            'cplsts3_stress',
            'cplsts4_stress',
            'cplsts6_stress',
            'cplsts8_stress',

            'cplstn3_strain',
            'cplstn4_strain',
            'cplstn6_strain',
            'cplstn8_strain',
            'cplsts3_strain',
            'cplsts4_strain',
            'cplsts6_strain',
            'cplsts8_strain',
        ]

        table_types += [
            # LAMA
            'eigenvalues',

            # HISADD
            #'convergence_history',

            # R1TABRG
            #'response1_table',

            # OEF - Forces - tCode=4 thermal=0
            'crod_force',
            'conrod_force',
            'ctube_force',

            # bar/beam/bend
            'cbend_force',

            'cbush_force',
            'coneax_force',
            'cdamp1_force',
            'cdamp2_force',
            'cdamp3_force',
            'cdamp4_force',
            'cgap_force',

            'cquad4_force',
            'cquad8_force',
            'cquadr_force',

            'ctria3_force',
            'ctria6_force',
            'ctriar_force',

            'cshear_force',
            #'cquad4_composite_force',
            #'cquad8_composite_force',
            #'cquadr_composite_force',
            #'ctria3_composite_force',
            #'ctria6_composite_force',
            #'ctriar_composite_force',

            'chexa_pressure_force',
            'cpenta_pressure_force',
            'ctetra_pressure_force',

            'celas1_force',
            'celas2_force',
            'celas3_force',
            'celas4_force',

            'cvisc_force',

            'force_VU',
            'force_VU_2D',

            #OEF - Fluxes - tCode=4 thermal=1
            'thermalLoad_CONV',

            #'thermalLoad_CHBDY',
            'chbdye_thermal_load',
            'chbdyg_thermal_load',
            'chbdyp_thermal_load',

            #'thermalLoad_1D',
            'crod_thermal_load',
            'cbeam_thermal_load',
            'ctube_thermal_load',
            'conrod_thermal_load',
            'cbar_thermal_load',
            'cbend_thermal_load',

            #'thermalLoad_2D_3D',
            'cquad4_thermal_load',
            'ctriax6_thermal_load',
            'cquad8_thermal_load',
            'ctria3_thermal_load',
            'ctria6_thermal_load',
            'ctetra_thermal_load',
            'cthexa_thermal_load',
            'cpenta_thermal_load',

            'thermalLoad_VU',
            'thermalLoad_VU_3D',
            'thermalLoad_VUBeam',
            #self.temperatureForces
        ]
        table_types += [
            # OES - CTRIAX6
            'ctriax_stress',
            'ctriax_strain',

            'cbush_stress',
            'cbush_strain',
            'cbush1d_stress_strain',

            # OES - nonlinear CROD/CONROD/CTUBE stress
            'nonlinear_crod_stress',
            'nonlinear_crod_strain',

            'nonlinear_ctube_stress',
            'nonlinear_ctube_strain',

            'nonlinear_conrod_stress',
            'nonlinear_conrod_strain',

            # OESNLXR - CTRIA3/CQUAD4 stress
            'nonlinear_cquad4_stress',
            'nonlinear_ctria3_stress',
            'nonlinear_cquad4_strain',
            'nonlinear_ctria3_strain',

            #'hyperelastic_plate_stress',
            'hyperelastic_cquad4_strain',

            # OES - CEALS1 224, CELAS3 225
            'nonlinear_celas1_stress',
            'nonlinear_celas3_stress',

            # OES - composite CTRIA3/CQUAD4 stress
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'cquadr_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',
            'ctriar_composite_stress',

            'cquad4_composite_strain',
            'cquad8_composite_strain',
            'cquadr_composite_strain',
            'ctria3_composite_strain',
            'ctria6_composite_strain',
            'ctriar_composite_strain',

            # OGS1 - grid point stresses
            'grid_point_stresses', # tCode=26
            'grid_point_volume_stresses',  # tCode=27

            # OEE - strain energy density
            #'strain_energy',  # tCode=18
            'cquad4_strain_energy',
            'cquad8_strain_energy',
            'cquadr_strain_energy',
            'cquadx_strain_energy',

            'ctria3_strain_energy',
            'ctria6_strain_energy',
            'ctriar_strain_energy',
            'ctriax_strain_energy',
            'ctriax6_strain_energy',

            'cshear_strain_energy',

            'ctetra_strain_energy',
            'cpenta_strain_energy',
            'chexa_strain_energy',
            'cpyram_strain_energy',

            'crod_strain_energy',
            'ctube_strain_energy',
            'conrod_strain_energy',

            'cbar_strain_energy',
            'cbeam_strain_energy',

            'cgap_strain_energy',
            'cbush_strain_energy',
            'celas1_strain_energy',
            'celas2_strain_energy',
            'celas3_strain_energy',
            'celas4_strain_energy',

            'cdum8_strain_energy',
            #'chexa8fd_strain_energy'
            'cbend_strain_energy',
            'dmig_strain_energy',
            'genel_strain_energy',


            # unused?
            'displacement_scaled_response_spectra_NRL',
            'displacement_scaled_response_spectra_ABS',
            'displacement_scaled_response_spectra_SRSS',
            'velocity_scaled_response_spectra_ABS',
            'acceleration_scaled_response_spectra_NRL',
            'acceleration_scaled_response_spectra_ABS',

        ]
        utables = unique(table_types)
        if len(table_types) != len(utables):
            msg = 'Non-unique tables: '
            for i, table_type in enumerate(table_types):
                if table_type in table_types[i+1:]:
                    msg += table_type + ', '
            raise AssertionError(msg)

        return table_types

    def _get_table_types_testing(self):
        """
        testing method...don't use
        """
        table_types = self.get_table_types()
        tables = object_attributes(self, 'public')
        tables = [table for table in tables
                  if isinstance(getattr(self, table), dict)
                  and table not in [
                      'card_count', 'data_code', 'element_mapper', 'iSubcaseNameMap',
                      'labels', 'subtitles', 'additional_matrices', 'matrices', 'subcase_key']]
        for table in tables:
            if self.make_geom:
                break
            assert table in table_types, table
        return table_types

    def get_f06_stats(self):
        return self.get_op2_stats()

    def get_op2_stats(self, short=False):
        """
        Gets info about the contents of the different attributes of the
        OP2 class.

        Example 1
        ---------
        >>> self.get_op2_stats()

        displacements[1]
          isubcase = 1
          type=RealDisplacementArray nnodes=72
          data: [t1, t2, t3, r1, r2, r3] shape=[1, 72, 6] dtype=float32
          gridTypes
          sort1
          lsdvmns = [1]

        spc_forces[1]
          isubcase = 1
          type=RealSPCForcesArray nnodes=72
          data: [t1, t2, t3, r1, r2, r3] shape=[1, 72, 6] dtype=float32
          gridTypes
          sort1
          lsdvmns = [1]

        ctetra_stress[1]
          type=RealSolidStressArray nelements=186 nnodes=930
          nodes_per_element=5 (including centroid)
          eType, cid
          data: [1, nnodes, 10] where 10=[oxx, oyy, ozz, txy, tyz, txz, o1, o2, o3, von_mises]
          data.shape = (1, 930, 10)
          element name: CTETRA
          sort1
          lsdvmns = [1]

        Example 1
        ---------
        >>> self.get_op2_stats(short=True)
        displacements[1]; RealDisplacementArray; [1, 72, 6]; [t1, t2, t3, r1, r2, r3]
        """
        def compare(key_value):
            key = key_value[0]
            if isinstance(key, (int, int32, int64, text_type, binary_type)):
                return key
            else:
                #print('key=%s type=%s' % (key, type(key)))
                #self.log.debug(type(key))
                return key[0]

        msg = []
        if self.grid_point_weight:
            msg += self.grid_point_weight.get_stats(short=short)

        table_types = self._get_table_types_testing()
        if short:
            no_data_classes = ['RealEigenvalues', 'ComplexEigenvalues', 'BucklingEigenvalues']
            for table_type in table_types:
                table = getattr(self, table_type)
                for isubcase, subcase in sorted(iteritems(table), key=compare):
                    class_name = subcase.__class__.__name__
                    if class_name in no_data_classes:
                        msg.append('%s[%r]\n' % (table_type, isubcase))
                    elif hasattr(subcase, 'data'):
                        #data = subcase.data
                        #shape = [int(i) for i in subcase.data.shape]
                        #headers = subcase.get_headers()
                        #headers_str = str(', '.join(headers))
                        #msg.append('%s[%s]; %s; %s; [%s]\n' % (
                        #table_type, isubcase, class_name, shape, headers_str))
                        msg.append('%s[%s]\n' % (table_type, isubcase))
                    elif hasattr(subcase, 'get_stats'):
                        msgi = '%s[%s] # unvectorized\n' % (table_type, isubcase)
                        msg.append(msgi)
                    else:
                        msgi = 'skipping %r %s[%s]\n' % (class_name, table_type, isubcase)
                        msg.append(msgi)
                        #raise RuntimeError(msgi)
        else:
            for table_type in table_types:
                table = getattr(self, table_type)
                try:
                    for isubcase, subcase in sorted(iteritems(table), key=compare):
                        class_name = subcase.__class__.__name__
                        if hasattr(subcase, 'get_stats'):
                            try:
                                stats = subcase.get_stats() # short=short
                            except:
                                msgi = 'errored reading %s %s[%s]\n\n' % (
                                    class_name, table_type, isubcase)
                                msg.append(msgi)
                                raise
                            else:
                                msg.append('%s[%s]\n' % (table_type, isubcase))
                                msg.extend(stats)
                                msg.append('\n')
                        else:
                            msgi = 'skipping %s %s[%s]\n\n' % (class_name, table_type, isubcase)
                            msg.append(msgi)
                            raise RuntimeError(msgi)
                except:
                    self.log.warning('type(table)=%s' % type(table))
                    self.log.warning(table)
                    raise

        for name, matrix in iteritems(self.matrices):
            msg.append('matrices[%s].shape = %s\n' % (name, matrix.data.shape))
        try:
            return ''.join(msg)
        except TypeError:
            for msgi in msg:
                print(msgi.rstrip())
                assert isinstance(msgi, string_types), msgi
        except UnicodeDecodeError:
            for msgi in msg:
                print(msgi.rstrip())
                assert isinstance(msgi, string_types), msgi
            raise

class Op2F06Attributes(OP2_F06_Common):
    def __init__(self):
        OP2_F06_Common.__init__(self)
