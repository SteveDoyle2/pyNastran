from __future__ import print_function
from collections import defaultdict
from six import string_types, binary_type, text_type
from numpy import unique, int32, int64

from pyNastran import is_release
from pyNastran.f06.f06_formatting import get_key0
from pyNastran.utils import object_attributes
from pyNastran.utils.numpy_utils import integer_types
from pyNastran.bdf.cards.base_card import deprecated
from pyNastran.bdf.case_control_deck import CaseControlDeck

from pyNastran.op2.result_objects.grid_point_weight import GridPointWeight
from pyNastran.op2.result_objects.design_response import Responses
from pyNastran.op2.result_objects.op2_results import Results


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
        self.params = {}
        self.table_names = []

        self.make_geom = False

        #: BDF Title
        self.title = None

        self.page_num = 1

        self.isubcases = []
        self.__objects_vector_init__()
        self.__objects_init__()
        self.__objects_common_init__()

    def has_result(self, result_name):
        """checks to see if a result exists"""
        if '.' in result_name:
            sline = result_name.split('.')
            if len(sline) != 2:
                msg = 'result_name=%r has too many dots; only 2 are allowed' % result_name
                raise RuntimeError(msg)

            #obj_names = sline[:-1]
            #storage_obj = self
            #for obj_name in obj_names:
                #print(obj_name)
                #storage_obj = getattr(storage_obj, obj_name)

            obj_name, result_name = result_name.split('.')
            try:
                storage_obj = getattr(self.op2_results, obj_name)
            except AttributeError:
                return False
            return hasattr(storage_obj, result_name)
        else:
            return hasattr(self, result_name)

    def get_result(self, result_name):
        """
        Getattr, but considers sub-objects

        Examples
        --------
        **Example 1**
        >>> self.eigenvectors = get_result('eigenvectors')

        **Example 1**
        >>> self.ato.displacements = get_result('ato.displacements')

        """
        if '.' in result_name:
            sline = result_name.split('.')
            if len(sline) != 2:
                raise RuntimeError('result_name=%r has too many dots; '
                                   'only 2 are allowed' % result_name)

            obj_name, result_name = result_name.split('.')

            storage_obj = getattr(self.op2_results, obj_name)
            storage_dict = getattr(storage_obj, result_name)
            return storage_dict
        else:
            try:
                storage_obj = getattr(self, result_name)
            except AttributeError:
                storage_obj = getattr(self.op2_results, result_name)
            return storage_obj

    def del_result(self, result_name):
        """
        delattr, but considers sub-objects
        """
        if '.' in result_name:
            sline = result_name.split('.')
            if len(sline) != 2:
                raise RuntimeError('result_name=%r has too many dots; '
                                   'only 2 are allowed' % result_name)
            obj_name, result_name = result_name.split('.')

            storage_obj = getattr(self.op2_results, obj_name)
            delattr(storage_obj, result_name)
        else:
            try:
                delattr(self, result_name)
            except AttributeError:
                storage_obj = delattr(self.op2_results, result_name)

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
        self.matrices = {}
        self.matdicts = {}
        #======================================================================
        # rods
        self.op2_results = Results()

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
        self.cbar_force_abs = {} # thermal=2
        self.cbar_force_nrl = {} # thermal=8

        self.cbar_stress = {}

        self.cbar_strain = {}

        self.cbar_force_10nodes = {}
        self.cbar_stress_10nodes = {}
        self.cbar_strain_10nodes = {}

        self.cbeam_force = {}
        self.cbeam_force_vu = {}

        self.cbeam_stress = {}
        self.cbeam_strain = {}

        #======================================================================
        self.cbend_stress = {}
        self.cbend_strain = {}
        self.cbend_force = {}

        #======================================================================
        # shells
        self.ctria3_force = {}
        self.ctria6_force = {}
        self.ctriar_force = {}
        self.cquad4_force = {}
        self.cquad8_force = {}
        self.cquadr_force = {}

        self.cquad4_composite_force_failure_indicies = {}
        self.cquad8_composite_force_failure_indicies = {}
        self.ctria3_composite_force_failure_indicies = {}
        self.ctria6_composite_force_failure_indicies = {}

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
        self.nonlinear_cbush1d_stress_strain = {}

        #======================================================================

    def __objects_common_init__(self):
        #: the date the job was run on
        self.date = None
        self.responses = Responses()

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
        self.displacements_scaled = {}    # tCode=1 thermal=8

        #: OUP

        self.displacement_scaled_response_spectra_nrl = {}  # thermal=8
        self.displacement_scaled_response_spectra_abs = {}  # thermal=2
        self.displacement_scaled_response_spectra_srss = {} # thermal=4
        #self.displacement_scaled_response_spectra_psd = {}
        #self.displacement_scaled_response_spectra_ato = {}
        #self.displacement_scaled_response_spectra_rms = {}
        #self.displacement_scaled_response_spectra_crm = {}
        #self.displacement_scaled_response_spectra_no = {}


        #: OUG - velocity
        self.velocities = {}              # tCode=10 thermal=0

        #self.velocity_scaled_response_spectra_nrl = {}
        self.velocity_scaled_response_spectra_abs = {}
        #self.velocity_scaled_response_spectra_psd = {}
        #self.velocity_scaled_response_spectra_ato = {}
        #self.velocity_scaled_response_spectra_rms = {}
        #self.velocity_scaled_response_spectra_crm = {}
        #self.velocity_scaled_response_spectra_no = {}

        #: OUG - acceleration
        self.accelerations = {}            # tCode=11 thermal=0

        self.acceleration_scaled_response_spectra_nrl = {}
        self.acceleration_scaled_response_spectra_abs = {}
        #self.acceleration_scaled_response_spectra_psd = {}
        #self.acceleration_scaled_response_spectra_ato = {}
        #self.acceleration_scaled_response_spectra_rms = {}
        #self.acceleration_scaled_response_spectra_crm = {}
        #self.acceleration_scaled_response_spectra_no = {}


        #: OUG - temperatures
        self.temperatures = {}           # tCode=1 thermal=1

        #: OUG - eigenvectors
        self.eigenvectors = {}            # tCode=7 thermal=0

        # OEF - Forces - tCode=4 thermal=0
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
        self.spc_forces = {}  # OQG1, tCode=3?
        self.spc_forces_v = {} # OQGV1
        self.spc_forces_scaled_response_spectra_nrl = {}

        self.mpc_forces = {}  # tCode=39
        self.mpc_forces_RAQCONS = {}
        self.mpc_forces_RAQEATC = {}

        # OQG - thermal forces
        self.thermal_gradient_and_flux = {}

        #: OGF - grid point forces
        self.grid_point_forces = {}  # tCode=19

        #: OGS1 - grid point stresses
        self.grid_point_surface_stresses = {}       # tCode=26
        self.grid_point_stresses_volume_direct = {}  # tCode=27
        self.grid_point_stresses_volume_principal = {} # tCode=28
        self.grid_point_stress_discontinuities = {}  # tCode=35

        #: OPG - summation of loads for each element
        self.load_vectors = {}       # OPG1; tCode=2  thermal=0
        self.load_vectors_v = {}     # OPGV!
        self.thermal_load_vectors = {}  # tCode=2  thermal=1
        self.applied_loads = {}       # tCode=19 thermal=0
        self.force_vectors = {}       # tCode=12 thermal=0


        #: OEE - strain energy density; tCode=18
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
        self.conm2_strain_energy = {}
        self.rbe1_strain_energy = {}
        self.rbe3_strain_energy = {}

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
                    msg = (
                        'bad compression check...\n'
                        'keys0=%s type(key0)=%s\n'
                        'res_key=%s type(res_key)=%s' % (
                            key0, type(key0), res_key, type(res_key))
                    )
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

    def remove_empty_results(self):
        table_names = self.get_table_types() + ['grid_point_weight', 'oload_resultant']
        for table_name in table_names:
            obj = self.get_result(table_name)
            if obj is None:
                obj = self.del_result(table_name)
            elif isinstance(obj, dict):
                if len(obj) == 0:
                    obj = self.del_result(table_name)
                    #delattr(self, table_name)
                else:
                    print('saving %s dict' % table_name)
            else:
                print('saving %s' % table_name)

    def get_table_types(self):
        """Gets the names of the results."""
        base = self.op2_results.get_table_types()

        table_types = base + [
            # OUG - displacement, temperatures, eigenvectors, velocity, acceleration
            'displacements', 'displacements_scaled',
            'temperatures',
            'eigenvectors',
            'velocities',
            'accelerations',

            # OQG - spc/mpc forces
            'spc_forces', 'spc_forces_v', 'spc_forces_scaled_response_spectra_nrl',
            'mpc_forces', 'mpc_forces_RAQCONS', 'mpc_forces_RAQEATC',
            'thermal_gradient_and_flux',

            # OGF - grid point forces
            'grid_point_forces',

            # OPG - summation of loads for each element '
            'load_vectors', 'load_vectors_v',
            'thermal_load_vectors',
            'applied_loads',
            'force_vectors',

            # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
            # OES - CELAS1/CELAS2/CELAS3/CELAS4 stress/strain
            'celas1_stress', 'celas2_stress', 'celas3_stress', 'celas4_stress',
            'celas1_strain', 'celas2_strain', 'celas3_strain', 'celas4_strain',

            # OES - isotropic CROD/CONROD/CTUBE stress/strain
            'crod_stress', 'conrod_stress', 'ctube_stress',
            'crod_strain', 'conrod_strain', 'ctube_strain',

            # OES - isotropic CBAR stress/strain
            'cbar_stress', 'cbar_strain',
            'cbar_force', 'cbar_force_abs', 'cbar_force_nrl',

            'cbar_stress_10nodes', 'cbar_strain_10nodes', 'cbar_force_10nodes',

            # OES - isotropic CBEAM stress/strain
            'cbeam_stress', 'cbeam_strain', 'cbeam_force', 'cbeam_force_vu',
            'nonlinear_cbeam_stress',
            #'nonlinear_cbeam_strain',

            # CBEND
            'cbend_stress', 'cbend_strain', 'cbend_force',

            # OES - isotropic CTRIA3/CQUAD4 stress
            'ctria3_stress', 'ctriar_stress', 'ctria6_stress',
            'cquadr_stress', 'cquad4_stress', 'cquad8_stress',

            # OES - isotropic CTRIA3/CQUAD4 strain
            'ctria3_strain', 'ctriar_strain', 'ctria6_strain',
            'cquadr_strain', 'cquad4_strain', 'cquad8_strain',

            # OES - isotropic CTETRA/CHEXA/CPENTA stress/strain
            'ctetra_stress', 'chexa_stress', 'cpenta_stress',
            'ctetra_strain', 'chexa_strain', 'cpenta_strain',

            # OES - CSHEAR stress/strain
            'cshear_stress', 'cshear_strain',

            # OES - GAPNL 86
            'nonlinear_cgap_stress',
            # OES - CBUSH 226
            'nonlinear_cbush_stress',

            'cplstn3_stress', 'cplstn4_stress', 'cplstn6_stress', 'cplstn8_stress',
            'cplsts3_stress', 'cplsts4_stress', 'cplsts6_stress', 'cplsts8_stress',

            'cplstn3_strain', 'cplstn4_strain', 'cplstn6_strain', 'cplstn8_strain',
            'cplsts3_strain', 'cplsts4_strain', 'cplsts6_strain', 'cplsts8_strain',
        ]

        table_types += [
            # PVT0
            'params',

            # LAMA
            'eigenvalues',

            # HISADD
            #'convergence_history',

            # R1TABRG
            #'response1_table',

            # OEF - Forces - tCode=4 thermal=0
            'crod_force', 'conrod_force', 'ctube_force',

            # bar/beam
            'cbush_force',
            'coneax_force',
            'cdamp1_force', 'cdamp2_force', 'cdamp3_force', 'cdamp4_force',
            'cgap_force',

            'cquad4_force', 'cquad8_force', 'cquadr_force',
            'ctria3_force', 'ctria6_force', 'ctriar_force',

            'cshear_force',
            #'cquad4_composite_force',
            #'cquad8_composite_force',
            #'cquadr_composite_force',
            #'ctria3_composite_force',
            #'ctria6_composite_force',
            #'ctriar_composite_force',

            'chexa_pressure_force', 'cpenta_pressure_force', 'ctetra_pressure_force',

            'celas1_force', 'celas2_force', 'celas3_force', 'celas4_force',
            'cvisc_force',

            'force_VU',
            #'force_VU_2D',
            'vu_quad_force', 'vu_tria_force',

            #OEF - Fluxes - tCode=4 thermal=1
            'conv_thermal_load',

            'cquad4_composite_force_failure_indicies',
            'cquad8_composite_force_failure_indicies',
            'ctria3_composite_force_failure_indicies',
            'ctria6_composite_force_failure_indicies',

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
            'ctriax_stress', 'ctriax_strain',

            'cbush_stress', 'cbush_strain',

            #'cbush_stress',
            #'cbush_strain',
            'cbush1d_stress_strain', 'nonlinear_cbush1d_stress_strain',

            # OES - nonlinear CROD/CONROD/CTUBE stress
            'nonlinear_crod_stress', 'nonlinear_crod_strain',

            'nonlinear_ctube_stress', 'nonlinear_ctube_strain',

            'nonlinear_conrod_stress', 'nonlinear_conrod_strain',

            # OESNLXR - CTRIA3/CQUAD4 stress
            'nonlinear_cquad4_stress', 'nonlinear_ctria3_stress',
            'nonlinear_cquad4_strain', 'nonlinear_ctria3_strain',

            #'hyperelastic_plate_stress',
            'hyperelastic_cquad4_strain',

            # OES - CEALS1 224, CELAS3 225
            'nonlinear_celas1_stress',
            'nonlinear_celas3_stress',

            # OES - composite CTRIA3/CQUAD4 stress
            'cquad4_composite_stress', 'cquad8_composite_stress', 'cquadr_composite_stress',
            'ctria3_composite_stress', 'ctria6_composite_stress', 'ctriar_composite_stress',

            'cquad4_composite_strain', 'cquad8_composite_strain', 'cquadr_composite_strain',
            'ctria3_composite_strain', 'ctria6_composite_strain', 'ctriar_composite_strain',

            # OGS1 - grid point stresses
            'grid_point_surface_stresses', # tCode=26
            'grid_point_stresses_volume_direct',  # tCode=27 # volume direct
            'grid_point_stresses_volume_principal', # tCode =28
            'grid_point_stress_discontinuities',  # tCode=35,

            # OEE - strain energy density
            #'strain_energy',  # tCode=18
            'cquad4_strain_energy', 'cquad8_strain_energy', 'cquadr_strain_energy',
            'cquadx_strain_energy',

            'ctria3_strain_energy', 'ctria6_strain_energy', 'ctriar_strain_energy',
            'ctriax_strain_energy', 'ctriax6_strain_energy',

            'cshear_strain_energy',

            'ctetra_strain_energy', 'cpenta_strain_energy',
            'chexa_strain_energy', 'cpyram_strain_energy',

            'crod_strain_energy', 'ctube_strain_energy', 'conrod_strain_energy',

            'cbar_strain_energy', 'cbeam_strain_energy',

            'cgap_strain_energy',
            'cbush_strain_energy',
            'celas1_strain_energy', 'celas2_strain_energy',
            'celas3_strain_energy', 'celas4_strain_energy',

            'cdum8_strain_energy',
            #'chexa8fd_strain_energy'
            'cbend_strain_energy',
            'dmig_strain_energy',
            'genel_strain_energy',
            'conm2_strain_energy',
            'rbe1_strain_energy', 'rbe3_strain_energy',

            # unused?
            'displacement_scaled_response_spectra_nrl',
            'displacement_scaled_response_spectra_abs',
            'displacement_scaled_response_spectra_srss',
            'velocity_scaled_response_spectra_abs',
            'acceleration_scaled_response_spectra_nrl',
            'acceleration_scaled_response_spectra_abs',
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
        skipped_attributes = [
            'card_count', 'data_code', 'element_mapper', 'isubcase_name_map',
            'labels', 'subtitles', 'additional_matrices', 'matrices', 'matdicts',
            'subcase_key', 'end_options', 'expected_times', 'generalized_tables',
            'op2_reader']

        table_types = self.get_table_types()
        tables = object_attributes(self, 'public')
        tables = [table for table in tables
                  if isinstance(getattr(self, table), dict)
                  and table not in skipped_attributes]
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
        msg += self.responses.get_stats(short=short)
        if self.grid_point_weight:
            msg += self.grid_point_weight.get_stats(short=short)

        table_types = self._get_table_types_testing()
        def _write_params(params):
            msg = []
            for key, param in sorted(params.items()):
                msg.append('PARAM[%s] = %s\n' % (key, param.values))
            return msg

        if short:
            no_data_classes = ['RealEigenvalues', 'ComplexEigenvalues', 'BucklingEigenvalues']
            for table_type in table_types:
                if table_type in ['params']:
                    msg.extend(_write_params(self.params))
                    continue
                elif table_type in ['gpdt', 'eqexin']:
                    obj = self.get_result(table_type)
                    if obj is None:
                        continue
                    stats = obj.get_stats(short=True)
                    msg.extend(stats)
                    continue

                table = self.get_result(table_type)
                for isubcase, subcase in sorted(table.items(), key=compare):
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
                    elif table_type == 'params':
                        msgi = str(subcase)
                    elif hasattr(subcase, 'get_stats'):
                        msgi = '%s[%s] # unvectorized\n' % (table_type, isubcase)
                        msg.append(msgi)
                    else:
                        msgi = 'skipping %r %s[%s]\n' % (class_name, table_type, isubcase)
                        msg.append(msgi)
                        #raise RuntimeError(msgi)
        else:
            for table_type in table_types:
                table = self.get_result(table_type)
                if table_type in ['params']:
                    msg.extend(_write_params(self.params))
                    continue
                elif table_type in ['gpdt', 'eqexin']:
                    obj = self.get_result(table_type)
                    if obj is None:
                        continue
                    stats = obj.get_stats(short=False)
                    msg.extend(stats)
                    continue

                try:
                    for isubcase, subcase in sorted(table.items(), key=compare):
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

        for unused_name, matrix in sorted(self.matrices.items()):
            #msg.append('matrices[%s].shape = %s\n' % (name, matrix.data.shape))
            msg.append(str(matrix) + '\n')

        for unused_name, matrix_dict in sorted(self.matdicts.items()):
            #msg.append('matrices[%s].shape = %s\n' % (name, matrix.data.shape))
            msg.append(str(matrix_dict) + '\n')
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
