from pyNastran.f06.tables.grid_point_weight import GridPointWeight
from collections import defaultdict
from numpy import unique


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
        self.accelerations = {}           # tCode=11 thermal=0

        # OEF - Forces - tCode=4 thermal=0
        # rods
        self.rodForces = {}

        self.barForces = {}
        self.bar100Forces = {}
        self.beamForces = {}
        self.bendForces = {}
        self.bushForces = {}
        self.coneAxForces = {}
        self.damperForces = {}
        self.gapForces = {}

        self.plateForces = {}
        self.plateForces2 = {}
        self.compositePlateForces = {}

        self.shearForces = {}

        self.solidPressureForces = {}
        self.springForces = {}
        self.viscForces = {}

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
        #: OES - CELAS1/CELAS2/CELAS3/CELAS4 stress
        self.celasStress = {}

        #: OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
        self.celasStrain = {}

        #: OES - CTRIAX6
        self.ctriaxStress = {}
        self.ctriaxStrain = {}

        #: OES - isotropic CROD/CONROD/CTUBE stress
        self.rodStress = {}

        #: OES - isotropic CROD/CONROD/CTUBE strain
        self.rodStrain = {}

        #: OES - nonlinear CROD/CONROD/CTUBE stress
        self.nonlinearRodStress = {}
        #: OES - nonlinear CROD/CONROD/CTUBE strain
        self.nonlinearRodStrain = {}
        #: OES - isotropic CBAR stress
        self.barStress = {}
        #: OES - isotropic CBAR strain
        self.barStrain = {}
        #: OES - isotropic CBEAM stress
        self.beamStress = {}
        #: OES - isotropic CBEAM strain
        self.beamStrain = {}
        #: OES - isotropic CBUSH stress
        self.bushStress = {}
        #: OES - isotropic CBUSH strain
        self.bushStrain = {}
         #: OES - isotropic CBUSH1D strain/strain
        self.bush1dStressStrain = {}

        #: OES - isotropic CTRIA3/CQUAD4 stress
        self.plateStress = {}
        #: OES - isotropic CTRIA3/CQUAD4 strain
        self.plateStrain = {}

        #: OESNLXR - CTRIA3/CQUAD4 stress
        self.nonlinearPlateStress = {}
        #: OESNLXR - CTRIA3/CQUAD4 strain
        self.nonlinearPlateStrain = {}
        self.hyperelasticPlateStress = {}
        self.hyperelasticPlateStrain = {}

        #: OES - isotropic CTETRA/CHEXA/CPENTA stress
        self.solidStress = {}

        #: OES - isotropic CTETRA/CHEXA/CPENTA strain
        self.solidStrain = {}

        #: OES - composite CTRIA3/CQUAD4 stress
        self.compositePlateStress = {}
        #: OES - composite CTRIA3/CQUAD4 strain
        self.compositePlateStrain = {}


        #: OES - CSHEAR stress
        self.shearStress = {}
        #: OES - CSHEAR strain
        self.shearStrain = {}
        #: OES - CELAS1 224, CELAS3 225,
        self.nonlinearSpringStress = {}
        #: OES - GAPNL 86
        self.nonlinearGapStress = {}
        #: OES - CBUSH 226
        self.nolinearBushStress = {}

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
                key0 = res_type.keys()[0]
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
            'celasStress',  # non-vectorized
            'celas1_stress',  # vectorized
            'celas2_stress',
            'celas3_stress',
            'celas4_stress',

            # OES - CELAS1/CELAS2/CELAS3/CELAS4 strain
            'celasStrain',  # non-vectorized
            'celas1_strain',  # vectorized
            'celas2_strain',
            'celas3_strain',
            'celas4_strain',

            # OES - isotropic CROD/CONROD/CTUBE stress
            'rodStress',  # non-vectorized
            'crod_stress',  # vectorized
            'conrod_stress',
            'ctube_stress',

            # OES - isotropic CROD/CONROD/CTUBE strain
            'rodStrain',  # non-vectorized
            'crod_strain',  # vectorized
            'conrod_strain',
            'ctube_strain',

            # OES - isotropic CBAR stress/strain
            'barStress',  # non-vectorized
            'barStrain',
            'cbar_stress',  # vectorized
            'cbar_strain',
            'cbar_force',

            # OES - isotropic CBEAM stress/strain
            'beamStress',  # non-vectorized
            'beamStrain',
            'cbeam_stress',  # vectorized
            'cbeam_strain',
            'cbeam_force',

            # OES - isotropic CBEAM strain

            # OES - isotropic CTRIA3/CQUAD4 stress
            'plateStress',  # non-vectorized
            'ctria3_stress',  # vectorized
            'ctriar_stress',
            'ctria6_stress',

            'cquadr_stress',
            'cquad4_stress',
            'cquad8_stress',

            # OES - isotropic CTRIA3/CQUAD4 strain
            'plateStrain',  # non-vectorized
            'ctria3_strain',  # vectorized
            'ctriar_strain',
            'ctria6_strain',

            'cquadr_strain',
            'cquad4_strain',
            'cquad8_strain',


            # OES - isotropic CTETRA/CHEXA/CPENTA stress
            'solidStress',  # non-vectorized
            'ctetra_stress',  # vectorized
            'chexa_stress',
            'cpenta_stress',

            # OES - isotropic CTETRA/CHEXA/CPENTA strain
            'solidStrain',  # non-vectorized
            'ctetra_strain',  # vectorized
            'chexa_strain',
            'cpenta_strain',

            # OES - CSHEAR stress
            'shearStress',
            'cshear_stress',
            # OES - CSHEAR strain
            'shearStrain',
            'cshear_strain',
            # OES - CEALS1 224, CELAS3 225
            'nonlinearSpringStress',
            # OES - GAPNL 86
            'nonlinearGapStress',
            # OES - CBUSH 226
            'nolinearBushStress',
        ]

        table_types += [
            # LAMA
            'eigenvalues',

            # OEF - Forces - tCode=4 thermal=0
            'rodForces',  # non-vectorized
            'crod_force',  # vectorized
            'conrod_force',
            'ctube_force',

            'barForces',
            'bar100Forces',
            'beamForces',
            'bendForces',
            'bushForces',
            'coneAxForces',
            'damperForces',
            'gapForces',

            'plateForces',
            'ctria3_force',
            'cquad4_force',

            'plateForces2',
            'shearForces',
            'solidPressureForces',
            'springForces',
            'viscForces',

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
            'ctriaxStress',
            'ctriaxStrain',

            'bushStress',
            'bushStrain',
            'bush1dStressStrain',

            # OES - nonlinear CROD/CONROD/CTUBE stress
            'nonlinearRodStress',
            'nonlinearRodStrain',

            # OESNLXR - CTRIA3/CQUAD4 stress
            'nonlinearPlateStress',
            'nonlinearPlateStrain',
            'hyperelasticPlateStress',
            'hyperelasticPlateStrain',

            # OES - composite CTRIA3/CQUAD4 stress
            'compositePlateStress',
            'cquad4_composite_stress',
            'cquad8_composite_stress',
            'ctria3_composite_stress',
            'ctria6_composite_stress',

            'compositePlateStrain',
            'cquad4_composite_strain',
            'cquad8_composite_strain',
            'ctria3_composite_strain',
            'ctria6_composite_strain',

            # OGS1 - grid point stresses
            'gridPointStresses',        # tCode=26
            'gridPointVolumeStresses',  # tCode=27

            # OEE - strain energy density
            'strainEnergy',  # tCode=18
        ]
        assert len(table_types) == len(unique(table_types))
        return table_types

    def get_f06_stats(self):
        return self.get_op2_stats()

    def get_op2_stats(self):
        """
        Gets info about the contents of the different attributes of the
        OP2 class.
        """
        table_types = self.get_table_types()
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

