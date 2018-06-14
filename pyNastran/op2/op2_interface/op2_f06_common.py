from __future__ import print_function
from collections import defaultdict
from six import iteritems, string_types, binary_type, text_type
from numpy import unique, int32, int64

from pyNastran import is_release
from pyNastran.op2.tables.grid_point_weight import GridPointWeight
from pyNastran.f06.f06_formatting import get_key0
from pyNastran.utils import object_attributes, integer_types
from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.bdf.case_control_deck import CaseControlDeck


class OP2_F06_Common(object):
    def __init__(self):
        #: a dictionary that maps an integer of the subcaseName to the
        #: subcase_id
        self.isubcase_name_map = {}
        self.generalized_tables = {}
        self.subtitles = defaultdict(list)
        self.case_control_deck = CaseControlDeck([], log=self.log)
        self.labels = {}
        self.expected_times = {}

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

        self.modal_contribution_crod_stress = {}
        self.modal_contribution_conrod_stress = {}
        self.modal_contribution_ctube_stress = {}
        self.modal_contribution_crod_strain = {}
        self.modal_contribution_conrod_strain = {}
        self.modal_contribution_ctube_strain = {}

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
        self.modal_contribution_celas1_stress = {}
        self.modal_contribution_celas2_stress = {}
        self.modal_contribution_celas3_stress = {}
        self.modal_contribution_celas4_stress = {}

        self.celas1_strain = {}
        self.celas2_strain = {}
        self.celas3_strain = {}
        self.celas4_strain = {}

        #======================================================================
        self.ctetra_stress = {}
        self.ctetra_stress_ato = {}
        self.ctetra_stress_crm = {}
        self.ctetra_stress_no = {}
        self.ctetra_stress_psd = {}
        self.ctetra_stress_rms = {}

        self.cpenta_stress = {}
        self.cpenta_stress_ato = {}
        self.cpenta_stress_crm = {}
        self.cpenta_stress_no = {}
        self.cpenta_stress_psd = {}
        self.cpenta_stress_rms = {}

        self.chexa_stress = {}
        self.chexa_stress_ato = {}
        self.chexa_stress_crm = {}
        self.chexa_stress_no = {}
        self.chexa_stress_psd = {}
        self.chexa_stress_rms = {}

        self.ctetra_strain = {}
        self.ctetra_strain_ato = {}
        self.ctetra_strain_crm = {}
        self.ctetra_strain_no = {}
        self.ctetra_strain_psd = {}
        self.ctetra_strain_rms = {}

        self.cpenta_strain = {}
        self.cpenta_strain_ato = {}
        self.cpenta_strain_crm = {}
        self.cpenta_strain_no = {}
        self.cpenta_strain_psd = {}
        self.cpenta_strain_rms = {}

        self.chexa_strain = {}
        self.chexa_strain_ato = {}
        self.chexa_strain_crm = {}
        self.chexa_strain_no = {}
        self.chexa_strain_psd = {}
        self.chexa_strain_rms = {}

        self.modal_contribution_ctetra_stress = {}
        self.modal_contribution_cpenta_stress = {}
        self.modal_contribution_chexa_stress = {}

        self.modal_contribution_ctetra_strain = {}
        self.modal_contribution_cpenta_strain = {}
        self.modal_contribution_chexa_strain = {}

        self.ctetra_stress_RASCONS = {}
        self.cpenta_stress_RASCONS = {}
        self.chexa_stress_RASCONS = {}

        self.ctetra_strain_RASCONS = {}
        self.cpenta_strain_RASCONS = {}
        self.chexa_strain_RASCONS = {}

        #======================================================================

        # bars/beams
        self.cbar_force = {}
        self.cbar_force_ato = {}
        self.cbar_force_crm = {}
        self.cbar_force_psd = {}
        self.cbar_force_rms = {}
        self.cbar_force_no = {}
        self.cbar_force_RAFCONS = {}
        self.cbar_force_RAFEATC = {}

        self.cbar_stress = {}
        self.cbar_stress_no = {}
        self.cbar_stress_ato = {}
        self.cbar_stress_crm = {}
        self.cbar_stress_psd = {}
        self.cbar_stress_rms = {}
        self.modal_contribution_cbar_stress = {}

        self.cbar_strain = {}
        self.cbar_strain_no = {}
        self.cbar_strain_ato = {}
        self.cbar_strain_crm = {}
        self.cbar_strain_psd = {}
        self.cbar_strain_rms = {}
        self.modal_contribution_cbar_strain = {}

        self.cbar_force_10nodes = {}
        self.cbar_stress_10nodes = {}
        self.cbar_strain_10nodes = {}

        self.cbeam_force = {}
        self.cbeam_force_ato = {}
        self.cbeam_force_crm = {}
        self.cbeam_force_psd = {}
        self.cbeam_force_rms = {}
        self.cbeam_force_no = {}

        self.cbeam_force_vu = {}

        self.cbeam_stress = {}
        self.cbeam_strain = {}
        self.modal_contribution_cbeam_stress = {}
        self.modal_contribution_cbeam_strain = {}

        #======================================================================
        # shells
        self.ctria3_force = {}
        self.ctria3_force_ato = {}
        self.ctria3_force_crm = {}
        self.ctria3_force_psd = {}
        self.ctria3_force_rms = {}
        self.ctria3_force_no = {}

        self.ctria6_force = {}
        self.ctria6_force_ato = {}
        self.ctria6_force_crm = {}
        self.ctria6_force_psd = {}
        self.ctria6_force_rms = {}
        self.ctria6_force_no = {}

        self.ctriar_force = {}
        self.ctriar_force_ato = {}
        self.ctriar_force_crm = {}
        self.ctriar_force_psd = {}
        self.ctriar_force_rms = {}
        self.ctriar_force_no = {}

        self.cquad4_force = {}
        self.cquad4_force_ato = {}
        self.cquad4_force_crm = {}
        self.cquad4_force_psd = {}
        self.cquad4_force_rms = {}
        self.cquad4_force_no = {}
        self.cquad4_force_RAFCONS = {}
        self.cquad4_force_RAFEATC = {}

        self.cquad8_force = {}
        self.cquad8_force_ato = {}
        self.cquad8_force_crm = {}
        self.cquad8_force_psd = {}
        self.cquad8_force_rms = {}
        self.cquad8_force_no = {}

        self.cquadr_force = {}
        self.cquadr_force_ato = {}
        self.cquadr_force_crm = {}
        self.cquadr_force_psd = {}
        self.cquadr_force_rms = {}
        self.cquadr_force_no = {}

        self.ctria3_stress = {}
        self.ctria6_stress = {}
        self.cquad4_stress = {}
        self.cquad8_stress = {}
        self.cquadr_stress = {}
        self.ctriar_stress = {}

        self.modal_contribution_ctria3_stress = {}
        self.modal_contribution_ctria6_stress = {}
        self.modal_contribution_cquad4_stress = {}
        self.modal_contribution_cquad8_stress = {}
        self.modal_contribution_cquadr_stress = {}
        self.modal_contribution_ctriar_stress = {}

        self.ctria3_stress_no = {}
        self.ctria6_stress_no = {}
        self.cquad4_stress_no = {}
        self.cquad8_stress_no = {}
        self.cquadr_stress_no = {}
        self.ctriar_stress_no = {}

        self.ctria3_stress_crm = {}
        self.ctria6_stress_crm = {}
        self.cquad4_stress_crm = {}
        self.cquad8_stress_crm = {}
        self.cquadr_stress_crm = {}
        self.ctriar_stress_crm = {}

        self.ctria3_stress_ato = {}
        self.ctria6_stress_ato = {}
        self.cquad4_stress_ato = {}
        self.cquad8_stress_ato = {}
        self.cquadr_stress_ato = {}
        self.ctriar_stress_ato = {}

        self.ctria3_stress_psd = {}
        self.ctria6_stress_psd = {}
        self.cquad4_stress_psd = {}
        self.cquad8_stress_psd = {}
        self.cquadr_stress_psd = {}
        self.ctriar_stress_psd = {}

        self.ctria3_stress_rms = {}
        self.ctria6_stress_rms = {}
        self.cquad4_stress_rms = {}
        self.cquad8_stress_rms = {}
        self.cquadr_stress_rms = {}
        self.ctriar_stress_rms = {}

        self.ctria3_stress_RASCONS = {}
        self.ctria6_stress_RASCONS = {}
        self.cquad4_stress_RASCONS = {}
        self.cquad8_stress_RASCONS = {}
        self.cquadr_stress_RASCONS = {}
        self.ctriar_stress_RASCONS = {}

        self.ctria3_strain = {}
        self.ctria6_strain = {}
        self.cquad4_strain = {}
        self.cquad8_strain = {}
        self.cquadr_strain = {}
        self.ctriar_strain = {}

        self.ctria3_strain_ato = {}
        self.ctria6_strain_ato = {}
        self.cquad4_strain_ato = {}
        self.cquad8_strain_ato = {}
        self.cquadr_strain_ato = {}
        self.ctriar_strain_ato = {}

        self.ctria3_strain_crm = {}
        self.ctria6_strain_crm = {}
        self.cquad4_strain_crm = {}
        self.cquad8_strain_crm = {}
        self.cquadr_strain_crm = {}
        self.ctriar_strain_crm = {}

        self.ctria3_strain_no = {}
        self.ctria6_strain_no = {}
        self.cquad4_strain_no = {}
        self.cquad8_strain_no = {}
        self.cquadr_strain_no = {}
        self.ctriar_strain_no = {}

        self.ctria3_strain_psd = {}
        self.ctria6_strain_psd = {}
        self.cquad4_strain_psd = {}
        self.cquad8_strain_psd = {}
        self.cquadr_strain_psd = {}
        self.ctriar_strain_psd = {}

        self.ctria3_strain_rms = {}
        self.ctria6_strain_rms = {}
        self.cquad4_strain_rms = {}
        self.cquad8_strain_rms = {}
        self.cquadr_strain_rms = {}
        self.ctriar_strain_rms = {}

        self.modal_contribution_ctria3_strain = {}
        self.modal_contribution_ctria6_strain = {}
        self.modal_contribution_cquad4_strain = {}
        self.modal_contribution_cquad8_strain = {}
        self.modal_contribution_cquadr_strain = {}
        self.modal_contribution_ctriar_strain = {}

        self.ctria3_strain_RASCONS = {}
        self.ctria6_strain_RASCONS = {}
        self.cquad4_strain_RASCONS = {}
        self.cquad8_strain_RASCONS = {}
        self.cquadr_strain_RASCONS = {}
        self.ctriar_strain_RASCONS = {}

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
        self.modal_contribution_cshear_stress = {}
        self.modal_contribution_cshear_strain = {}
        self.modal_contribution_cshear_force = {}

        #: OES - CBEAM 94
        self.nonlinear_cbeam_stress = {}

        # bushing
        self.cbush_stress = {}
        self.cbush_strain = {}
        self.nonlinear_cbush_stress = {}  # CBUSH 226
        self.modal_contribution_cbush_stress = {}
        self.modal_contribution_cbush_strain = {}

        self.cbush1d_stress_strain = {}
        self.nonlinear_cbush1d_stress_strain = {}

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

        #: self.frequencies already exists as a BDF object
        #: but we need this for the FOL frequencies for the MONPNT1 and MONPNT3
        self._frequencies = None

        #: ESE
        self.eigenvalues = {}

        #: OUG - displacement
        self.displacements = {}           # tCode=1 thermal=0
        self.displacements_PSD = {}        # random
        self.displacements_ATO = {}        # random
        self.displacements_RMS = {}        # random
        self.displacements_CRM = {}        # random
        self.displacements_NO = {}         # random
        self.displacements_scaled = {}    # tCode=1 thermal=8
        self.displacements_ROUGV1 = {}

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
        self.velocities_PSD = {}
        self.velocities_ATO = {}
        self.velocities_RMS = {}
        self.velocities_CRM = {}
        self.velocities_NO = {}
        self.velocities_ROUGV1 = {}

        #self.velocity_scaled_response_spectra_NRL = {}
        self.velocity_scaled_response_spectra_ABS = {}
        #self.velocity_scaled_response_spectra_PSD = {}
        #self.velocity_scaled_response_spectra_ATO = {}
        #self.velocity_scaled_response_spectra_RMS = {}
        #self.velocity_scaled_response_spectra_CRM = {}
        #self.velocity_scaled_response_spectra_NO = {}

        #: OUG - acceleration
        self.accelerations = {}            # tCode=11 thermal=0
        self.accelerations_PSD = {}
        self.accelerations_ATO = {}
        self.accelerations_RMS = {}
        self.accelerations_CRM = {}
        self.accelerations_NO = {}
        self.accelerations_ROUGV1 = {}

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
        self.eigenvectors_ROUGV1 = {}

        # OEF - Forces - tCode=4 thermal=0

        self.cbend_force = {}

        self.cbush_force = {}
        self.cbush_force_ato = {}
        self.cbush_force_psd = {}
        self.cbush_force_crm = {}
        self.cbush_force_rms = {}
        self.cbush_force_no = {}
        self.cbush_force_RAFCONS = {}
        self.cbush_force_RAFEATC = {}

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

        #self.force_VU_2D = {}
        self.vu_quad_force = {}
        self.vu_tria_force = {}

        #OEF - Fluxes - tCode=4 thermal=1
        self.conv_thermal_load = {}

        #self.thermalLoad_CHBDY = {}
        self.chbdye_thermal_load = {}
        self.chbdyg_thermal_load = {}
        self.chbdyp_thermal_load = {}
        self.chbdye_thermal_load_flux = {}
        self.chbdyg_thermal_load_flux = {}
        self.chbdyp_thermal_load_flux = {}

        #self.thermalLoad_1D = {}
        self.crod_thermal_load = {}
        self.cbeam_thermal_load = {}
        self.ctube_thermal_load = {}
        self.conrod_thermal_load = {}
        self.cbar_thermal_load = {}
        self.cbend_thermal_load = {}
        self.crod_thermal_load_flux = {}
        self.cbeam_thermal_load_flux = {}
        self.ctube_thermal_load_flux = {}
        self.conrod_thermal_load_flux = {}
        self.cbar_thermal_load_flux = {}
        self.cbend_thermal_load_flux = {}

        #self.thermalLoad_2D_3D = {}
        self.cquad4_thermal_load = {}
        self.ctriax6_thermal_load = {}
        self.cquad8_thermal_load = {}
        self.ctria3_thermal_load = {}
        self.ctria6_thermal_load = {}
        self.ctetra_thermal_load = {}
        self.chexa_thermal_load = {}
        self.cpenta_thermal_load = {}

        self.cquad4_thermal_load_flux = {}
        self.ctriax6_thermal_load_flux = {}
        self.cquad8_thermal_load_flux = {}
        self.ctria3_thermal_load_flux = {}
        self.ctria6_thermal_load_flux = {}
        self.ctetra_thermal_load_flux = {}
        self.chexa_thermal_load_flux = {}
        self.cpenta_thermal_load_flux = {}

        self.thermalLoad_VU = {}
        self.thermalLoad_VU_3D = {}
        self.vu_beam_thermal_load = {}
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
        self.spc_forces_PSD = {}
        self.spc_forces_ATO = {}
        self.spc_forces_RMS = {}
        self.spc_forces_CRM = {}
        self.spc_forces_NO = {}

        self.mpc_forces = {}  # tCode=39
        self.mpc_forces_PSD = {}
        self.mpc_forces_ATO = {}
        self.mpc_forces_RMS = {}
        self.mpc_forces_CRM = {}
        self.mpc_forces_NO = {}
        self.mpc_forces_RAQCONS = {}
        self.mpc_forces_RAQEATC = {}

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
        self.load_vectors_ATO = {}
        self.load_vectors_CRM = {}
        self.load_vectors_NO = {}
        self.load_vectors_PSD = {}
        self.load_vectors_RMS = {}


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
        """Gets the names of the results."""
        table_types = [
            # OUG - displacement
            'displacements',
            'displacements_PSD',
            'displacements_ATO',
            'displacements_RMS',
            'displacements_CRM',
            'displacements_NO',
            'displacements_scaled',
            'displacements_ROUGV1',

            # OUG - temperatures
            'temperatures',

            # OUG - eigenvectors
            'eigenvectors',
            'eigenvectors_RADCONS',
            'eigenvectors_RADEFFM',
            'eigenvectors_RADEATC',
            'eigenvectors_ROUGV1',


            # OUG - velocity
            'velocities',
            'velocities_PSD',
            'velocities_ATO',
            'velocities_RMS',
            'velocities_CRM',
            'velocities_NO',
            'velocities_ROUGV1',

            # OUG - acceleration
            'accelerations',
            'accelerations_PSD',
            'accelerations_ATO',
            'accelerations_RMS',
            'accelerations_CRM',
            'accelerations_NO',
            'accelerations_ROUGV1',

            # OQG - spc/mpc forces
            'spc_forces', 'spc_forces_PSD', 'spc_forces_ATO', 'spc_forces_RMS', 'spc_forces_CRM', 'spc_forces_NO',
            'spc_forces_scaled_response_spectra_NRL',

            'mpc_forces', 'mpc_forces_PSD', 'mpc_forces_ATO', 'mpc_forces_RMS', 'mpc_forces_CRM', 'mpc_forces_NO',
            'mpc_forces_RAQCONS', 'mpc_forces_RAQEATC',

            'thermal_gradient_and_flux',

            # OGF - grid point forces
            'grid_point_forces',

            # OPG - summation of loads for each element '
            'load_vectors', 'load_vectors_ATO', 'load_vectors_CRM', 'load_vectors_NO', 'load_vectors_PSD', 'load_vectors_RMS',
            'thermal_load_vectors',
            'applied_loads',
            'force_vectors',

            # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
            # OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
            'celas1_stress', 'modal_contribution_celas1_stress',
            'celas2_stress', 'modal_contribution_celas2_stress',
            'celas3_stress', 'modal_contribution_celas3_stress',
            'celas4_stress', 'modal_contribution_celas4_stress',

            # OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
            'celas1_strain',
            'celas2_strain',
            'celas3_strain',
            'celas4_strain',

            # OES - isotropic CROD/CONROD/CTUBE stress
            'crod_stress', 'modal_contribution_crod_stress',
            'conrod_stress', 'modal_contribution_conrod_stress',
            'ctube_stress', 'modal_contribution_ctube_stress',

            # OES - isotropic CROD/CONROD/CTUBE strain
            'crod_strain', 'modal_contribution_crod_strain',
            'conrod_strain', 'modal_contribution_conrod_strain',
            'ctube_strain', 'modal_contribution_ctube_strain',

            # OES - isotropic CBAR stress/strain
            'cbar_stress', 'cbar_stress_ato', 'cbar_stress_crm', 'cbar_stress_no', 'cbar_stress_psd', 'cbar_stress_rms', 'modal_contribution_cbar_stress',
            'cbar_strain', 'cbar_strain_ato', 'cbar_strain_crm', 'cbar_strain_no', 'cbar_strain_psd', 'cbar_strain_rms', 'modal_contribution_cbar_strain',
            'cbar_force', 'cbar_force_ato', 'cbar_force_crm', 'cbar_force_psd', 'cbar_force_rms', 'cbar_force_no',
            'cbar_force_RAFCONS', 'cbar_force_RAFEATC',

            'cbar_stress_10nodes',
            'cbar_strain_10nodes',
            'cbar_force_10nodes',

            # OES - isotropic CBEAM stress/strain
            'cbeam_stress', 'modal_contribution_cbeam_stress',
            'cbeam_strain', 'modal_contribution_cbeam_strain',
            'cbeam_force', 'cbeam_force_ato', 'cbeam_force_crm', 'cbeam_force_psd', 'cbeam_force_rms', 'cbeam_force_no',
            'cbeam_force_vu',
            'nonlinear_cbeam_stress',
            #'nonlinear_cbeam_strain',


            # OES - isotropic CTRIA3/CQUAD4 stress
            'ctria3_stress', 'ctria3_stress_ato', 'ctria3_stress_crm', 'ctria3_stress_no', 'ctria3_stress_psd', 'ctria3_stress_rms', 'modal_contribution_ctria3_stress', 'ctria3_stress_RASCONS',
            'ctriar_stress', 'ctriar_stress_ato', 'ctriar_stress_crm', 'ctriar_stress_no', 'ctriar_stress_psd', 'ctriar_stress_rms', 'modal_contribution_ctriar_stress', 'ctriar_stress_RASCONS',
            'ctria6_stress', 'ctria6_stress_ato', 'ctria6_stress_crm', 'ctria6_stress_no', 'ctria6_stress_psd', 'ctria6_stress_rms', 'modal_contribution_ctria6_stress', 'ctria6_stress_RASCONS',

            'cquadr_stress', 'cquadr_stress_ato', 'cquadr_stress_crm', 'cquadr_stress_no', 'cquadr_stress_psd', 'cquadr_stress_rms', 'modal_contribution_cquadr_stress', 'cquadr_stress_RASCONS',
            'cquad4_stress', 'cquad4_stress_ato', 'cquad4_stress_crm', 'cquad4_stress_no', 'cquad4_stress_psd', 'cquad4_stress_rms', 'modal_contribution_cquad4_stress', 'cquad4_stress_RASCONS',
            'cquad8_stress', 'cquad8_stress_ato', 'cquad8_stress_crm', 'cquad8_stress_no', 'cquad8_stress_psd', 'cquad8_stress_rms', 'modal_contribution_cquad8_stress', 'cquad8_stress_RASCONS',

            # OES - isotropic CTRIA3/CQUAD4 strain
            'ctria3_strain', 'ctria3_strain_ato', 'ctria3_strain_crm', 'ctria3_strain_no', 'ctria3_strain_psd', 'ctria3_strain_rms', 'modal_contribution_ctria3_strain', 'ctria3_strain_RASCONS',
            'ctriar_strain', 'ctriar_strain_ato', 'ctriar_strain_crm', 'ctriar_strain_no', 'ctriar_strain_psd', 'ctriar_strain_rms', 'modal_contribution_ctriar_strain', 'ctriar_strain_RASCONS',
            'ctria6_strain', 'ctria6_strain_ato', 'ctria6_strain_crm', 'ctria6_strain_no', 'ctria6_strain_psd', 'ctria6_strain_rms', 'modal_contribution_ctria6_strain', 'ctria6_strain_RASCONS',

            'cquadr_strain', 'cquadr_strain_ato', 'cquadr_strain_crm', 'cquadr_strain_no', 'cquadr_strain_psd', 'cquadr_strain_rms', 'modal_contribution_cquadr_strain', 'cquadr_strain_RASCONS',
            'cquad4_strain', 'cquad4_strain_ato', 'cquad4_strain_crm', 'cquad4_strain_no', 'cquad4_strain_psd', 'cquad4_strain_rms', 'modal_contribution_cquad4_strain', 'cquad4_strain_RASCONS',
            'cquad8_strain', 'cquad8_strain_ato', 'cquad8_strain_crm', 'cquad8_strain_no', 'cquad8_strain_psd', 'cquad8_strain_rms', 'modal_contribution_cquad8_strain', 'cquad8_strain_RASCONS',


            # OES - isotropic CTETRA/CHEXA/CPENTA stress
            'ctetra_stress', 'ctetra_stress_ato', 'ctetra_stress_crm', 'ctetra_stress_no', 'ctetra_stress_psd', 'ctetra_stress_rms', 'modal_contribution_ctetra_stress', 'ctetra_stress_RASCONS',
            'chexa_stress',  'chexa_stress_ato',  'chexa_stress_crm',  'chexa_stress_no',  'chexa_stress_psd',  'chexa_stress_rms',  'modal_contribution_chexa_stress', 'chexa_stress_RASCONS',
            'cpenta_stress', 'cpenta_stress_ato', 'cpenta_stress_crm', 'cpenta_stress_no', 'cpenta_stress_psd', 'cpenta_stress_rms', 'modal_contribution_cpenta_stress', 'cpenta_stress_RASCONS',

            # OES - isotropic CTETRA/CHEXA/CPENTA strain
            'ctetra_strain', 'ctetra_strain_ato', 'ctetra_strain_crm', 'ctetra_strain_no', 'ctetra_strain_psd', 'ctetra_strain_rms', 'modal_contribution_ctetra_strain', 'ctetra_strain_RASCONS',
            'chexa_strain',  'chexa_strain_ato',  'chexa_strain_crm',  'chexa_strain_no',  'chexa_strain_psd',  'chexa_strain_rms',  'modal_contribution_chexa_strain', 'chexa_strain_RASCONS',
            'cpenta_strain', 'cpenta_strain_ato', 'cpenta_strain_crm', 'cpenta_strain_no', 'cpenta_strain_psd', 'cpenta_strain_rms', 'modal_contribution_cpenta_strain', 'cpenta_strain_RASCONS',

            # OES - CSHEAR stress/strain
            'cshear_stress', 'modal_contribution_cshear_stress',
            'cshear_strain', 'modal_contribution_cshear_strain',

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

            'cbush_force', 'cbush_force_ato', 'cbush_force_crm', 'cbush_force_psd', 'cbush_force_rms', 'cbush_force_no',
            'cbush_force_RAFCONS', 'cbush_force_RAFEATC',
            'coneax_force',
            'cdamp1_force',
            'cdamp2_force',
            'cdamp3_force',
            'cdamp4_force',
            'cgap_force',

            'cquad4_force', 'cquad4_force_ato', 'cquad4_force_crm', 'cquad4_force_psd', 'cquad4_force_rms', 'cquad4_force_no',
            'cquad8_force', 'cquad8_force_ato', 'cquad8_force_crm', 'cquad8_force_psd', 'cquad8_force_rms', 'cquad8_force_no',
            'cquadr_force', 'cquadr_force_ato', 'cquadr_force_crm', 'cquadr_force_psd', 'cquadr_force_rms', 'cquadr_force_no',

            'ctria3_force', 'ctria3_force_ato', 'ctria3_force_crm', 'ctria3_force_psd', 'ctria3_force_rms', 'ctria3_force_no',
            'ctria6_force', 'ctria6_force_ato', 'ctria6_force_crm', 'ctria6_force_psd', 'ctria6_force_rms', 'ctria6_force_no',
            'ctriar_force', 'ctriar_force_ato', 'ctriar_force_crm', 'ctriar_force_psd', 'ctriar_force_rms', 'ctriar_force_no',
            'cquad4_force_RAFCONS', 'cquad4_force_RAFEATC',

            'cshear_force', 'modal_contribution_cshear_force',
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
            #'force_VU_2D',
            'vu_quad_force',
            'vu_tria_force',

            #OEF - Fluxes - tCode=4 thermal=1
            'conv_thermal_load',

            #'thermalLoad_CHBDY',
            'chbdye_thermal_load', 'chbdye_thermal_load_flux',
            'chbdyg_thermal_load', 'chbdyg_thermal_load_flux',
            'chbdyp_thermal_load', 'chbdyp_thermal_load_flux',

            #'thermalLoad_1D',
            'crod_thermal_load', 'crod_thermal_load_flux',
            'cbeam_thermal_load', 'cbeam_thermal_load_flux',
            'ctube_thermal_load', 'ctube_thermal_load_flux',
            'conrod_thermal_load', 'conrod_thermal_load_flux',
            'cbar_thermal_load', 'cbar_thermal_load_flux',
            'cbend_thermal_load', 'cbend_thermal_load_flux',

            #'thermalLoad_2D_3D',
            'cquad4_thermal_load', 'cquad4_thermal_load_flux',
            'ctriax6_thermal_load', 'ctriax6_thermal_load_flux',
            'cquad8_thermal_load', 'cquad8_thermal_load_flux',
            'ctria3_thermal_load', 'ctria3_thermal_load_flux',
            'ctria6_thermal_load', 'ctria6_thermal_load_flux',
            'ctetra_thermal_load', 'ctetra_thermal_load_flux',
            'chexa_thermal_load', 'chexa_thermal_load_flux',
            'cpenta_thermal_load', 'cpenta_thermal_load_flux',

            'thermalLoad_VU',
            'thermalLoad_VU_3D',
            'vu_beam_thermal_load',
            #self.temperatureForces
        ]
        table_types += [
            # OES - CTRIAX6
            'ctriax_stress',
            'ctriax_strain',

            'cbush_stress', 'modal_contribution_cbush_stress',
            'cbush_strain', 'modal_contribution_cbush_strain',
            'cbush1d_stress_strain', 'nonlinear_cbush1d_stress_strain',

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
        """testing method...don't use"""
        table_types = self.get_table_types()
        tables = object_attributes(self, 'public')
        tables = [table for table in tables
                  if isinstance(getattr(self, table), dict)
                  and table not in [
                      'card_count', 'data_code', 'element_mapper', 'isubcase_name_map',
                      'labels', 'subtitles', 'additional_matrices', 'matrices', 'subcase_key',
                      'end_options', 'expected_times', 'generalized_tables']]
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

        Examples
        --------
        ***Detailed OP2 Stats***
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

        ***Appreviated OP2 Stats***
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

        for name, matrix in sorted(iteritems(self.matrices)):
            #msg.append('matrices[%s].shape = %s\n' % (name, matrix.data.shape))
            msg.append(str(matrix) + '\n')
        try:
            return ''.join(msg)
        except TypeError:
            for msgi in msg:
                print('TypeError...%r' % msgi.rstrip())
                assert isinstance(msgi, string_types), msgi
        except UnicodeDecodeError:
            for msgi in msg:
                print('UnicodeDecodeError...%r' % msgi.rstrip())
                assert isinstance(msgi, string_types), msgi
            raise

class Op2F06Attributes(OP2_F06_Common):
    def __init__(self):
        OP2_F06_Common.__init__(self)
