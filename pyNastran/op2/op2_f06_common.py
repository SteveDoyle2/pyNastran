from six import iteritems
from collections import defaultdict
from numpy import unique

from pyNastran.f06.tables.grid_point_weight import GridPointWeight
from pyNastran.f06.f06_formatting import get_key0
from pyNastran.utils import object_attributes


class OP2_F06_Common(object):
    def __init__(self):
        #: a dictionary that maps an integer of the subcaseName to the
        #: subcaseID
        self.iSubcaseNameMap = {}
        self.subtitles = defaultdict(list)
        self.labels = {}

        #: BDF Title
        self.Title = None

        self.page_num = 1

        self.iSubcases = []
        self.__objects_vector_init__()
        self.__objects_init__()
        self.__objects_common_init__()


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

        self.cbeam_force = {}
        self.cbeam_stress = {}
        self.cbeam_strain = {}

        #======================================================================
        # shells
        self.ctria3_force = {}
        self.cquad4_force = {}

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

        #: Grid Point Weight Table
        #: create with:
        #:   PARAM   GRDPNT    0  (required for F06/OP2)
        #:   PARAM   POSTEXT YES  (required for OP2)
        self.grid_point_weight = GridPointWeight()
        self.oload_resultant = None

        #: ESE
        self.eigenvalues = {}

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
        self.scaledDisplacements = {}     # tCode=1 thermal=8

        #: OUP

        self.displacement_scaled_response_spectra_NRL = {}
        self.displacement_scaled_response_spectra_ABS = {}
        #self.displacement_scaled_response_spectra_PSD = {}
        #self.displacement_scaled_response_spectra_ATO = {}
        #self.displacement_scaled_response_spectra_RMS = {}
        #self.displacement_scaled_response_spectra_CRM = {}
        #self.displacement_scaled_response_spectra_NO = {}

        #self.velocity_scaled_response_spectra_NRL = {}
        self.velocity_scaled_response_spectra_ABS = {}
        #self.velocity_scaled_response_spectra_PSD = {}
        #self.velocity_scaled_response_spectra_ATO = {}
        #self.velocity_scaled_response_spectra_RMS = {}
        #self.velocity_scaled_response_spectra_CRM = {}
        #self.velocity_scaled_response_spectra_NO = {}

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

        #: OUG - velocity
        self.velocities = {}              # tCode=10 thermal=0

        #: OUG - acceleration
        self.accelerations = {}            # tCode=11 thermal=0

        # OEF - Forces - tCode=4 thermal=0

        self.cbar100_force = {}
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

        self.plateForces = {}
        self.plateForces2 = {}
        self.compositePlateForces = {}


        #self.solidPressureForces = {}
        self.chexa_pressure_force = {}
        self.cpenta_pressure_force = {}
        self.ctetra_pressure_force = {}

        self.cvisc_force = {}

        self.force_VU = {}
        self.force_VU_2D = {}

        #OEF - Fluxes - tCode=4 thermal=1
        self.thermalLoad_CONV = {}
        self.thermalLoad_CHBDY = {}
        self.thermalLoad_1D = {}
        self.thermalLoad_2D_3D = {}
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

         #: OES - isotropic CBUSH1D strain/strain
        #self.bush1dStressStrain = {}

        #: OESNLXR - CTRIA3/CQUAD4 stress
        #self.nonlinearPlateStress = {}
        #: OESNLXR - CTRIA3/CQUAD4 strain
        #self.nonlinearPlateStrain = {}
        self.hyperelasticPlateStress = {}
        self.hyperelasticPlateStrain = {}

        #: OES - composite CTRIA3/CQUAD4 stress
        #self.compositePlateStress = {}
        self.nonlinear_cquad4_stress = {}
        self.nonlinear_ctria3_stress = {}
        #: OES - composite CTRIA3/CQUAD4 strain
        #self.compositePlateStrain = {}
        self.nonlinear_cquad4_strain = {}
        self.nonlinear_ctria3_strain = {}


        #: OES - CELAS1 224, CELAS3 225,
        self.nonlinearSpringStress = {}
        #: OES - GAPNL 86
        #self.nonlinearGapStress = {}
        self.nonlinear_cgap_stress = {}

        # OQG - spc/mpc forces
        self.spcForces = {}  # tCode=3?
        self.mpcForces = {}  # tCode=39

        # OQG - thermal forces
        self.thermalGradientAndFlux = {}

        #: OGF - grid point forces
        self.gridPointForces = {}  # tCode=19

        #: OGS1 - grid point stresses
        self.gridPointStresses = {}       # tCode=26
        self.gridPointVolumeStresses = {}  # tCode=27

        #: OPG - summation of loads for each element
        self.loadVectors = {}       # tCode=2  thermal=0
        self.thermalLoadVectors = {}  # tCode=2  thermal=1
        self.appliedLoads = {}       # tCode=19 thermal=0
        self.forceVectors = {}       # tCode=12 thermal=0

        #: OEE - strain energy density
        self.strainEnergy = {}  # tCode=18

    def _get_result_length(self, res_types, res_key):
        res_length = 0
        for res_type in res_types:
            if res_key in res_type:
                result = res_type[res_key]
                res_length = max(len(result.__class__.__name__), res_length)
                continue
            elif len(res_type) != 0:
                key0 = get_key0(res_type)
                result = res_type[key0]
                class_name = result.__class__.__name__
                print('%s - results not found...key=%s' % (class_name, res_key))
            else:  # empty result
                pass
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
            'scaledDisplacements',

            # OUG - temperatures
            'temperatures',

            # OUG - eigenvectors
            'eigenvectors',

            # OUG - velocity
            'velocities',

            # OUG - acceleration
            'accelerations',

            # OQG - spc/mpc forces
            'spcForces',
            'mpcForces',
            'thermalGradientAndFlux',

            # OGF - grid point forces
            'gridPointForces',

            # OPG - summation of loads for each element
            'loadVectors',
            'thermalLoadVectors',
            'appliedLoads',
            'forceVectors',

            # OES - tCode=5 thermal=0 s_code=0,1 (stress/strain)
            # OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
            #'celasStress',  # non-vectorized
            'celas1_stress',  # vectorized
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
            #'rodStrain',  # non-vectorized
            'crod_strain',
            'conrod_strain',
            'ctube_strain',

            # OES - isotropic CBAR stress/strain
            'cbar_stress',
            'cbar_strain',
            'cbar_force',
            'cbar100_force',

            # OES - isotropic CBEAM stress/strain
            'cbeam_stress',
            'cbeam_strain',
            'cbeam_force',
            'nonlinear_cbeam_stress',
            #'nonlinear_cbeam_strain',


            # OES - isotropic CTRIA3/CQUAD4 stress
            #'plateStress',  # non-vectorized
            'ctria3_stress',
            'ctriar_stress',
            'ctria6_stress',

            'cquadr_stress',
            'cquad4_stress',
            'cquad8_stress',

            # OES - isotropic CTRIA3/CQUAD4 strain
            #'plateStrain',  # non-vectorized
            'ctria3_strain',  # vectorized
            'ctriar_strain',
            'ctria6_strain',

            'cquadr_strain',
            'cquad4_strain',
            'cquad8_strain',


            # OES - isotropic CTETRA/CHEXA/CPENTA stress
            #'solidStress',  # non-vectorized
            'ctetra_stress',  # vectorized
            'chexa_stress',
            'cpenta_stress',

            # OES - isotropic CTETRA/CHEXA/CPENTA strain
            #'solidStrain',  # non-vectorized
            'ctetra_strain',  # vectorized
            'chexa_strain',
            'cpenta_strain',

            # OES - CSHEAR stress/strain
            'cshear_stress',
            'cshear_strain',
            # OES - CEALS1 224, CELAS3 225
            'nonlinearSpringStress',
            # OES - GAPNL 86
            'nonlinear_cgap_stress',
            # OES - CBUSH 226
            'nonlinear_cbush_stress',
        ]

        table_types += [
            # LAMA
            'eigenvalues',

            # OEF - Forces - tCode=4 thermal=0
            'crod_force',
            'conrod_force',
            'ctube_force',

            # bar/beam/bend
            #'bar100Force',
            'cbend_force',

            'cbush_force',
            'coneax_force',
            'cdamp1_force',
            'cdamp2_force',
            'cdamp3_force',
            'cdamp4_force',
            'cgap_force',

            'plateForces',
            'ctria3_force',
            'cquad4_force',

            'plateForces2',
            #'shearForces',
            'cshear_force',
            'compositePlateForces',

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
            'thermalLoad_CHBDY',
            'thermalLoad_1D',
            'thermalLoad_2D_3D',
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
            #'nonlinearPlateStress',
            #'nonlinearPlateStrain',
            'nonlinear_cquad4_stress',
            'nonlinear_ctria3_stress',
            'nonlinear_cquad4_strain',
            'nonlinear_ctria3_strain',

            'hyperelasticPlateStress',
            'hyperelasticPlateStrain',

            # OES - composite CTRIA3/CQUAD4 stress
            #'compositePlateStress',
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'cquadr_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',
            'ctriar_composite_stress',

            #'compositePlateStrain',
            'cquad4_composite_strain',
            'cquad8_composite_strain',
            'cquadr_composite_strain',
            'ctria3_composite_strain',
            'ctria6_composite_strain',
            'ctriar_composite_strain',

            # OGS1 - grid point stresses
            'gridPointStresses',        # tCode=26
            'gridPointVolumeStresses',  # tCode=27

            # OEE - strain energy density
            'strainEnergy',  # tCode=18

            # unused?
            'displacement_scaled_response_spectra_NRL',
            'displacement_scaled_response_spectra_ABS',
            'velocity_scaled_response_spectra_ABS',
            'acceleration_scaled_response_spectra_NRL',
            'acceleration_scaled_response_spectra_ABS',

        ]
        assert len(table_types) == len(unique(table_types))
        return table_types

    def _get_table_types_testing(self):
        """
        testing method...don't use
        """
        table_types = self.get_table_types()
        tables = object_attributes(self, 'public')
        tables = [table for table in tables
                  if isinstance(getattr(self, table), dict)
                  and table not in ['card_count', 'data_code', 'element_mapper', 'iSubcaseNameMap',
                  'labels', 'subtitles']]
        for table in tables:
            assert table in table_types, table
        return table_types

    def get_f06_stats(self):
        return self.get_op2_stats()

    def get_op2_stats(self):
        """
        Gets info about the contents of the different attributes of the
        OP2 class.
        """
        table_types = self._get_table_types_testing()
        msg = []
        for table_type in table_types:
            table = getattr(self, table_type)
            for isubcase, subcase in sorted(iteritems(table)):
                if hasattr(subcase, 'get_stats'):
                    msg.append('op2.%s[%s]\n' % (table_type, isubcase))
                    msg.extend(subcase.get_stats())
                    msg.append('\n')
                else:
                    msg.append('skipping %s op2.%s[%s]\n\n' % (subcase.__class__.__name__, table_type, isubcase))
                    #raise RuntimeError('skipping %s op2.%s[%s]\n\n' % (subcase.__class__.__name__, table_type, isubcase))
        return ''.join(msg)

