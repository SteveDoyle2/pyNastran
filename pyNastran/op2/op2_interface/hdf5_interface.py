# pylint: disable=C0301
"""
defines:
 model = load_op2_from_h5(h5_filename, log=None)
 export_op2_to_hdf5(hdf5_filename, op2_model)

 model = load_op2_from_hdf5(hdf5_filename, combine=True, log=None)
 model = load_op2_from_hdf5_file(model, h5_file, log, debug=False)
 export_op2_to_hdf5_file(hdf5_filename, op2_model)
 export_op2_to_hdf5_file(hdf5_file, op2_model)

"""
from __future__ import (nested_scopes, generators, division, absolute_import,
                        print_function, unicode_literals)
import os
from six import b

from six import iteritems, PY3, binary_type
import numpy as np
import h5py

import pyNastran
from pyNastran.op2.op2 import OP2

from pyNastran.op2.tables.lama_eigenvalues.lama_objects import RealEigenvalues, ComplexEigenvalues, BucklingEigenvalues
from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray, ComplexDisplacementArray
from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray, ComplexVelocityArray
from pyNastran.op2.tables.oug.oug_accelerations import RealAccelerationArray, ComplexAccelerationArray
from pyNastran.op2.tables.oug.oug_eigenvectors import RealEigenvectorArray, ComplexEigenvectorArray
from pyNastran.op2.tables.oug.oug_temperatures import RealTemperatureArray
from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import RealThermalVelocityVectorArray


from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import (
    RealLoadVectorArray, ComplexLoadVectorArray,
    RealTemperatureVectorArray,
)
from pyNastran.op2.tables.oqg_constraintForces.oqg_thermal_gradient_and_flux import (
    #RealTemperatureGradientAndFlux,
    RealTemperatureGradientAndFluxArray)

from pyNastran.op2.tables.opg_appliedLoads.opnl_force_vector import RealForceVectorArray#, ComplexForceVectorArray

from pyNastran.op2.tables.oqg_constraintForces.oqg_spc_forces import RealSPCForcesArray, ComplexSPCForcesArray
from pyNastran.op2.tables.oqg_constraintForces.oqg_mpc_forces import RealMPCForcesArray, ComplexMPCForcesArray
#from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import RealLoadVectorArray, ComplexLoadVectorArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_plates import RealPlateStressArray, RealPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_rods import RealRodStressArray, RealRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars import RealBarStressArray, RealBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bars100 import RealBar10NodesStressArray, RealBar10NodesStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_beams import RealBeamStressArray, RealBeamStrainArray, RealNonlinearBeamStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_shear import RealShearStressArray, RealShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids import RealSolidStressArray, RealSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import RealSpringStressArray, RealSpringStrainArray, RealNonlinearSpringStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray, RealCompositePlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush import RealBushStressArray, RealBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_gap import NonlinearGapStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_triax import RealTriaxStressArray #, RealTriaxStrainArray

from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_rod import RealNonlinearRodArray
from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import RealNonlinearPlateArray
#from pyNastran.op2.tables.oes_stressStrain.oes_hyperelastic import HyperelasticQuadArray

from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray, ComplexBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_beams import ComplexBeamStressArray, ComplexBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush import ComplexCBushStressArray, ComplexCBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush1d import ComplexCBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import ComplexPlateStressArray, ComplexPlateStrainArray, ComplexTriaxStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray, ComplexRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import ComplexShearStressArray, ComplexShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray, ComplexSpringStrainArray

from pyNastran.op2.tables.oes_stressStrain.random.oes_rods import RandomRodStressArray, RandomRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bars import RandomBarStressArray, RandomBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_beams import RandomBeamStressArray, RandomBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates import RandomPlateStressArray, RandomPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids import RandomSolidStressArray, RandomSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_shear import RandomShearStressArray, RandomShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_composite_plates import RandomCompositePlateStressArray, RandomCompositePlateStrainArray


from pyNastran.op2.tables.ogs_grid_point_stresses.ogs_surface_stresses import GridPointStressesArray, GridPointStressesVolumeArray


from pyNastran.op2.tables.oee_energy.oee_objects import RealStrainEnergyArray, ComplexStrainEnergyArray
from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray, ComplexGridPointForcesArray
from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    RealBendForceArray, RealCBar100ForceArray, RealCBarForceArray, RealCBeamForceArray,
    RealCBeamForceVUArray, RealCBushForceArray, RealCGapForceArray, RealConeAxForceArray,
    RealCShearForceArray, RealDamperForceArray, RealForceVU2DArray, RealPlateBilinearForceArray,
    RealPlateForceArray, RealRodForceArray, RealSolidPressureForceArray, RealSpringForceArray,
    RealViscForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexCBarForceArray, ComplexCBeamForceArray, ComplexCBeamForceVUArray,
    ComplexCBendForceArray, ComplexCBushForceArray, ComplexCShearForceArray,
    ComplexDamperForceArray, ComplexPlate2ForceArray, ComplexPlateForceArray,
    ComplexRodForceArray, ComplexSolidPressureForceArray, ComplexSpringForceArray,
    ComplexViscForceArray,
)
from pyNastran.op2.tables.oef_forces.oef_thermal_objects import (
    RealChbdyHeatFluxArray, # RealConvHeatFluxArray,
    Real1DHeatFluxArray,
    #RealElementTableArray, RealHeatFluxVUArray, RealHeatFluxVUBeamArray,
    #RealHeatFluxVU3DArray,
    HeatFlux_2D_3DArray,
)
#from pyNastran.op2.tables.oqg_constraintForces.oqg_thermal_gradient_and_flux import RealTemperatureGradientAndFluxArray
from pyNastran.utils import print_bad_path

def _cast(h5_result_attr):
    """converts the h5py type back into the OP2 type"""
    if h5_result_attr is None:
        return None

    if len(h5_result_attr.shape) == 0:
        return np.array(h5_result_attr).tolist()
        #raise NotImplementedError(h5_result_attr.dtype)
    return np.array(h5_result_attr)

TABLE_OBJ_MAP = {
    'displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'velocities' : (RealVelocityArray, ComplexVelocityArray, RealThermalVelocityVectorArray),
    'accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'eigenvectors' : (RealEigenvectorArray, ComplexEigenvectorArray),
    'spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'mpc_forces' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'load_vectors' : (RealLoadVectorArray, ComplexLoadVectorArray),
    'force_vectors' : (RealForceVectorArray, ),

    'displacement_scaled_response_spectra_ABS' : (RealDisplacementArray, ComplexDisplacementArray),
    'velocity_scaled_response_spectra_ABS' : (RealVelocityArray, ComplexVelocityArray),
    'spc_forces_scaled_response_spectra_ABS' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'acceleration_scaled_response_spectra_ABS' : (RealAccelerationArray, ComplexAccelerationArray),

    'displacement_scaled_response_spectra_NRL' : (RealDisplacementArray, ComplexDisplacementArray),
    'spc_forces_scaled_response_spectra_NRL' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'acceleration_scaled_response_spectra_NRL' : (RealAccelerationArray, ComplexAccelerationArray),

    'displacement_scaled_response_spectra_SRSS' : (RealDisplacementArray, ComplexDisplacementArray),
    'acceleration_scaled_response_spectra_SRSS' : (RealAccelerationArray, ComplexAccelerationArray),

    'displacements_NO' : (RealDisplacementArray, ComplexDisplacementArray),
    'displacements_ATO' : (RealDisplacementArray, ComplexDisplacementArray),
    'displacements_CRM' : (RealDisplacementArray, ComplexDisplacementArray),
    'displacements_PSD' : (RealDisplacementArray, ComplexDisplacementArray),
    'displacements_RMS' : (RealDisplacementArray, ComplexDisplacementArray),

    'velocities_NO' : (RealVelocityArray, ComplexVelocityArray),
    'velocities_ATO' : (RealVelocityArray, ComplexVelocityArray),
    'velocities_CRM' : (RealVelocityArray, ComplexVelocityArray),
    'velocities_PSD' : (RealVelocityArray, ComplexVelocityArray),
    'velocities_RMS' : (RealVelocityArray, ComplexVelocityArray),

    'accelerations_NO' : (RealAccelerationArray, ComplexAccelerationArray),
    'accelerations_ATO' : (RealAccelerationArray, ComplexAccelerationArray),
    'accelerations_CRM' : (RealAccelerationArray, ComplexAccelerationArray),
    'accelerations_PSD' : (RealAccelerationArray, ComplexAccelerationArray),
    'accelerations_RMS' : (RealAccelerationArray, ComplexAccelerationArray),

    'spc_forces_NO' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'spc_forces_ATO' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'spc_forces_CRM' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'spc_forces_PSD' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'spc_forces_RMS' : (RealSPCForcesArray, ComplexSPCForcesArray),

    'mpc_forces_NO' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'mpc_forces_ATO' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'mpc_forces_CRM' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'mpc_forces_PSD' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'mpc_forces_RMS' : (RealMPCForcesArray, ComplexMPCForcesArray),

    'celas1_stress' : (RealSpringStressArray, ComplexSpringStressArray),
    'celas1_stress_ato' : (RealSpringStressArray, ),
    'celas1_stress_crm' : (RealSpringStressArray, ),
    'celas1_stress_psd' : (RealSpringStressArray, ),
    'celas1_stress_rms' : (RealSpringStressArray, ),
    'celas1_stress_no' : (RealSpringStressArray, ),

    'celas2_stress' : (RealSpringStressArray, ComplexSpringStressArray),
    'celas2_stress_ato' : (RealSpringStressArray, ),
    'celas2_stress_crm' : (RealSpringStressArray, ),
    'celas2_stress_psd' : (RealSpringStressArray, ),
    'celas2_stress_rms' : (RealSpringStressArray, ),
    'celas2_stress_no' : (RealSpringStressArray, ),

    'celas3_stress' : (RealSpringStressArray, ComplexSpringStressArray),
    'celas3_stress_ato' : (RealSpringStressArray, ),
    'celas3_stress_crm' : (RealSpringStressArray, ),
    'celas3_stress_psd' : (RealSpringStressArray, ),
    'celas3_stress_rms' : (RealSpringStressArray, ),
    'celas3_stress_no' : (RealSpringStressArray, ),

    'celas4_stress' : (RealSpringStressArray, ComplexSpringStressArray),
    'celas4_stress_ato' : (RealSpringStressArray, ),
    'celas4_stress_crm' : (RealSpringStressArray, ),
    'celas4_stress_psd' : (RealSpringStressArray, ),
    'celas4_stress_rms' : (RealSpringStressArray, ),
    'celas4_stress_no' : (RealSpringStressArray, ),

    'celas1_strain' : (RealSpringStrainArray, ComplexSpringStrainArray),
    'celas1_strain_ato' : (RealSpringStrainArray, ),
    'celas1_strain_crm' : (RealSpringStrainArray, ),
    'celas1_strain_psd' : (RealSpringStrainArray, ),
    'celas1_strain_rms' : (RealSpringStrainArray, ),
    'celas1_strain_no' : (RealSpringStrainArray, ),

    'celas2_strain' : (RealSpringStrainArray, ComplexSpringStrainArray),
    'celas2_strain_ato' : (RealSpringStrainArray, ),
    'celas2_strain_crm' : (RealSpringStrainArray, ),
    'celas2_strain_psd' : (RealSpringStrainArray, ),
    'celas2_strain_rms' : (RealSpringStrainArray, ),
    'celas2_strain_no' : (RealSpringStrainArray, ),

    'celas3_strain' : (RealSpringStrainArray, ComplexSpringStrainArray),
    'celas3_strain_ato' : (RealSpringStrainArray, ),
    'celas3_strain_crm' : (RealSpringStrainArray, ),
    'celas3_strain_psd' : (RealSpringStrainArray, ),
    'celas3_strain_rms' : (RealSpringStrainArray, ),
    'celas3_strain_no' : (RealSpringStrainArray, ),

    'celas4_strain' : (RealSpringStrainArray, ComplexSpringStrainArray),
    'celas4_strain_ato' : (RealSpringStrainArray, ),
    'celas4_strain_crm' : (RealSpringStrainArray, ),
    'celas4_strain_psd' : (RealSpringStrainArray, ),
    'celas4_strain_rms' : (RealSpringStrainArray, ),
    'celas4_strain_no' : (RealSpringStrainArray, ),

    'celas1_force' : (RealSpringForceArray, ComplexSpringForceArray),
    'celas2_force' : (RealSpringForceArray, ComplexSpringForceArray),
    'celas3_force' : (RealSpringForceArray, ComplexSpringForceArray),
    'celas4_force' : (RealSpringForceArray, ComplexSpringForceArray),

    'cdamp1_force' : (RealDamperForceArray, ComplexDamperForceArray),
    'cdamp2_force' : (RealDamperForceArray, ComplexDamperForceArray),
    'cdamp3_force' : (RealDamperForceArray, ComplexDamperForceArray),
    'cdamp4_force' : (RealDamperForceArray, ComplexDamperForceArray),

    'cvisc_force' : (RealViscForceArray, ComplexViscForceArray),

    'crod_stress' : (RealRodStressArray, ComplexRodStressArray),
    'crod_stress_ato' : (RandomRodStressArray, ),
    'crod_stress_crm' : (RandomRodStressArray, ),
    'crod_stress_psd' : (RandomRodStressArray, ),
    'crod_stress_rms' : (RandomRodStressArray, ),
    'crod_stress_no' : (RandomRodStressArray, ),

    'conrod_stress' : (RealRodStressArray, ComplexRodStressArray),
    'conrod_stress_ato' : (RandomRodStressArray, ),
    'conrod_stress_crm' : (RandomRodStressArray, ),
    'conrod_stress_psd' : (RandomRodStressArray, ),
    'conrod_stress_rms' : (RandomRodStressArray, ),
    'conrod_stress_no' : (RandomRodStressArray, ),

    'ctube_stress' : (RealRodStressArray, ComplexRodStressArray),
    'ctube_stress_ato' : (RandomRodStressArray, ),
    'ctube_stress_crm' : (RandomRodStressArray, ),
    'ctube_stress_psd' : (RandomRodStressArray, ),
    'ctube_stress_rms' : (RandomRodStressArray, ),
    'ctube_stress_no' : (RandomRodStressArray, ),

    'crod_strain' : (RealRodStrainArray, ComplexRodStrainArray),
    'crod_strain_ato' : (RandomRodStrainArray, ),
    'crod_strain_crm' : (RandomRodStrainArray, ),
    'crod_strain_psd' : (RandomRodStrainArray, ),
    'crod_strain_rms' : (RandomRodStrainArray, ),
    'crod_strain_no' : (RandomRodStrainArray, ),

    'conrod_strain' : (RealRodStrainArray, ComplexRodStrainArray),
    'conrod_strain_ato' : (RandomRodStrainArray, ),
    'conrod_strain_crm' : (RandomRodStrainArray, ),
    'conrod_strain_psd' : (RandomRodStrainArray, ),
    'conrod_strain_rms' : (RandomRodStrainArray, ),
    'conrod_strain_no' : (RandomRodStrainArray, ),

    'ctube_strain' : (RealRodStrainArray, ComplexRodStrainArray),
    'ctube_strain_ato' : (RandomRodStrainArray, ),
    'ctube_strain_crm' : (RandomRodStrainArray, ),
    'ctube_strain_psd' : (RandomRodStrainArray, ),
    'ctube_strain_rms' : (RandomRodStrainArray, ),
    'ctube_strain_no' : (RandomRodStrainArray, ),

    'crod_force' : (RealRodForceArray, ComplexRodForceArray),
    'conrod_force' : (RealRodForceArray, ComplexRodForceArray),
    'ctube_force' : (RealRodForceArray, ComplexRodForceArray),

    'cbar_stress' : (RealBarStressArray, ComplexBarStressArray),
    'cbar_stress_ato' : (RandomBarStressArray, ),
    'cbar_stress_crm' : (RandomBarStressArray, ),
    'cbar_stress_psd' : (RandomBarStressArray, ),
    'cbar_stress_no' : (RandomBarStressArray, ),
    'cbar_stress_rms' : (RandomBarStressArray, ),

    'cbar_strain' : (RealBarStrainArray, ComplexBarStrainArray),
    'cbar_strain_ato' : (RandomBarStrainArray, ),
    'cbar_strain_crm' : (RandomBarStrainArray, ),
    'cbar_strain_psd' : (RandomBarStrainArray, ),
    'cbar_strain_no' : (RandomBarStrainArray, ),
    'cbar_strain_rms' : (RandomBarStrainArray, ),

    'cbar_force' : (RealCBarForceArray, ComplexCBarForceArray),
    'cbar_force_ato' : (RealCBarForceArray, ),
    'cbar_force_crm' : (RealCBarForceArray, ),
    'cbar_force_psd' : (RealCBarForceArray, ),
    'cbar_force_no' : (RealCBarForceArray, ),
    'cbar_force_rms' : (RealCBarForceArray, ),

    'cbar_force_10nodes' : (RealCBar100ForceArray, ),
    'cbar_stress_10nodes' : (RealBar10NodesStressArray, ),
    'cbar_strain_10nodes' : (RealBar10NodesStrainArray, ),

    'cbeam_stress' : (RealBeamStressArray, ComplexBeamStressArray),
    'cbeam_stress_ato' : (RandomBeamStressArray, ),
    'cbeam_stress_crm' : (RandomBeamStressArray, ),
    'cbeam_stress_psd' : (RandomBeamStressArray, ),
    'cbeam_stress_rms' : (RandomBeamStressArray, ),
    'cbeam_stress_no' : (RandomBeamStressArray, ),

    'cbeam_strain' : (RealBeamStrainArray, ComplexBeamStrainArray),
    'cbeam_strain_ato' : (RandomBeamStrainArray, ),
    'cbeam_strain_crm' : (RandomBeamStrainArray, ),
    'cbeam_strain_psd' : (RandomBeamStrainArray, ),
    'cbeam_strain_rms' : (RandomBeamStrainArray, ),
    'cbeam_strain_no' : (RandomBeamStrainArray, ),

    'cbeam_force' : (RealCBeamForceArray, ComplexCBeamForceArray),

    'cquad4_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'cquad4_stress_ato' : (RandomPlateStressArray, ),
    'cquad4_stress_crm' : (RandomPlateStressArray, ),
    'cquad4_stress_psd' : (RandomPlateStressArray, ),
    'cquad4_stress_no' : (RandomPlateStressArray, ),
    'cquad4_stress_rms' : (RandomPlateStressArray, ),

    'ctria3_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'cquad4_stress_ato' : (RandomPlateStressArray, ),
    'cquad4_stress_crm' : (RandomPlateStressArray, ),
    'cquad4_stress_psd' : (RandomPlateStressArray, ),
    'cquad4_stress_no' : (RandomPlateStressArray, ),
    'cquad4_stress_rms' : (RandomPlateStressArray, ),

    'cquad8_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'cquad8_stress_ato' : (RandomPlateStressArray, ),
    'cquad8_stress_crm' : (RandomPlateStressArray, ),
    'cquad8_stress_psd' : (RandomPlateStressArray, ),
    'cquad8_stress_no' : (RandomPlateStressArray, ),
    'cquad8_stress_rms' : (RandomPlateStressArray, ),

    'ctria6_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'ctria6_stress_ato' : (RandomPlateStressArray, ),
    'ctria6_stress_crm' : (RandomPlateStressArray, ),
    'ctria6_stress_psd' : (RandomPlateStressArray, ),
    'ctria6_stress_no' : (RandomPlateStressArray, ),
    'ctria6_stress_rms' : (RandomPlateStressArray, ),

    'ctriar_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'ctriar_stress_ato' : (RandomPlateStressArray, ),
    'ctriar_stress_crm' : (RandomPlateStressArray, ),
    'ctriar_stress_psd' : (RandomPlateStressArray, ),
    'ctriar_stress_no' : (RandomPlateStressArray, ),
    'ctriar_stress_rms' : (RandomPlateStressArray, ),

    'cquadr_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'cquadr_stress_ato' : (RandomPlateStressArray, ),
    'cquadr_stress_crm' : (RandomPlateStressArray, ),
    'cquadr_stress_psd' : (RandomPlateStressArray, ),
    'cquadr_stress_no' : (RandomPlateStressArray, ),
    'cquadr_stress_rms' : (RandomPlateStressArray, ),

    'cquad4_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'cquad4_strain_ato' : (RandomPlateStrainArray, ),
    'cquad4_strain_crm' : (RandomPlateStrainArray, ),
    'cquad4_strain_psd' : (RandomPlateStrainArray, ),
    'cquad4_strain_no' : (RandomPlateStrainArray, ),
    'cquad4_strain_rms' : (RandomPlateStrainArray, ),

    'ctria3_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ctria3_strain_ato' : (RandomPlateStrainArray, ),
    'ctria3_strain_crm' : (RandomPlateStrainArray, ),
    'ctria3_strain_psd' : (RandomPlateStrainArray, ),
    'ctria3_strain_no' : (RandomPlateStrainArray, ),
    'ctria3_strain_rms' : (RandomPlateStrainArray, ),

    'cquad8_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'cquad8_strain_ato' : (RandomPlateStrainArray, ),
    'cquad8_strain_crm' : (RandomPlateStrainArray, ),
    'cquad8_strain_psd' : (RandomPlateStrainArray, ),
    'cquad8_strain_no' : (RandomPlateStrainArray, ),
    'cquad8_strain_rms' : (RandomPlateStrainArray, ),

    'ctria6_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ctria6_strain_ato' : (RandomPlateStrainArray, ),
    'ctria6_strain_crm' : (RandomPlateStrainArray, ),
    'ctria6_strain_psd' : (RandomPlateStrainArray, ),
    'ctria6_strain_no' : (RandomPlateStrainArray, ),
    'ctria6_strain_rms' : (RandomPlateStrainArray, ),

    'ctriar_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ctriar_strain_ato' : (RandomPlateStrainArray, ),
    'ctriar_strain_crm' : (RandomPlateStrainArray, ),
    'ctriar_strain_psd' : (RandomPlateStrainArray, ),
    'ctriar_strain_no' : (RandomPlateStrainArray, ),
    'ctriar_strain_rms' : (RandomPlateStrainArray, ),

    'cquadr_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'cquadr_strain_ato' : (RandomPlateStrainArray, ),
    'cquadr_strain_crm' : (RandomPlateStrainArray, ),
    'cquadr_strain_psd' : (RandomPlateStrainArray, ),
    'cquadr_strain_no' : (RandomPlateStrainArray, ),
    'cquadr_strain_rms' : (RandomPlateStrainArray, ),

    'ctria3_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'ctria3_stress_ato' : (RandomPlateStressArray, ),
    'ctria3_stress_crm' : (RandomPlateStressArray, ),
    'ctria3_stress_psd' : (RandomPlateStressArray, ),
    'ctria3_stress_no' : (RandomPlateStressArray, ),
    'ctria3_stress_rms' : (RandomPlateStressArray, ),

    'ctria3_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ctria3_strain_ato' : (RandomPlateStrainArray, ),
    'ctria3_strain_crm' : (RandomPlateStrainArray, ),
    'ctria3_strain_psd' : (RandomPlateStrainArray, ),
    'ctria3_strain_no' : (RandomPlateStrainArray, ),
    'ctria3_strain_rms' : (RandomPlateStrainArray, ),

    'cquad43_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'cquad4_stress_ato' : (RandomPlateStressArray, ),
    'cquad4_stress_crm' : (RandomPlateStressArray, ),
    'cquad4_stress_psd' : (RandomPlateStressArray, ),
    'cquad4_stress_no' : (RandomPlateStressArray, ),
    'cquad4_stress_rms' : (RandomPlateStressArray, ),

    'cquad4_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'cquad4_strain_ato' : (RandomPlateStrainArray, ),
    'cquad4_strain_crm' : (RandomPlateStrainArray, ),
    'cquad4_strain_psd' : (RandomPlateStrainArray, ),
    'cquad4_strain_no' : (RandomPlateStrainArray, ),
    'cquad4_strain_rms' : (RandomPlateStrainArray, ),

    'ctriax_stress' : (RealTriaxStressArray, ComplexTriaxStressArray,),

    'cquad4_composite_stress' : (RealCompositePlateStressArray, None),
    'ctria3_composite_stress' : (RealCompositePlateStressArray, None),
    'ctria6_composite_stress' : (RealCompositePlateStressArray, None),
    'cquad8_composite_stress' : (RealCompositePlateStressArray, None),

    'cquad4_composite_strain' : (RealCompositePlateStrainArray, None),
    'ctria3_composite_strain' : (RealCompositePlateStrainArray, None),
    'ctria6_composite_strain' : (RealCompositePlateStrainArray, None),
    'cquad8_composite_strain' : (RealCompositePlateStrainArray, None),

    'cshear_stress' : (RealShearStressArray, ComplexShearStressArray),
    'cshear_stress_ato' : (RandomShearStressArray, ),
    'cshear_stress_crm' : (RandomShearStressArray, ),
    'cshear_stress_psd' : (RandomShearStressArray, ),
    'cshear_stress_no' : (RandomShearStressArray, ),
    'cshear_stress_rms' : (RandomShearStressArray, ),

    'cshear_strain' : (RealShearStrainArray, ComplexShearStrainArray),
    'cshear_strain_ato' : (RandomShearStrainArray, ),
    'cshear_strain_crm' : (RandomShearStrainArray, ),
    'cshear_strain_psd' : (RandomShearStrainArray, ),
    'cshear_strain_no' : (RandomShearStrainArray, ),
    'cshear_strain_rms' : (RandomShearStrainArray, ),

    'cshear_force' : (RealCShearForceArray, ComplexCShearForceArray),

    'coneax_force' : (RealConeAxForceArray,),

    'ctetra_stress' : (RealSolidStressArray, ComplexSolidStressArray),
    'ctetra_stress_ato' : (RandomSolidStressArray, ),
    'ctetra_stress_crm' : (RandomSolidStressArray, ),
    'ctetra_stress_psd' : (RandomSolidStressArray, ),
    'ctetra_stress_no' : (RandomSolidStressArray, ),
    'ctetra_stress_rms' : (RandomSolidStressArray, ),

    'cpenta_stress' : (RealSolidStressArray, ComplexSolidStressArray),
    'cpenta_stress_ato' : (RandomSolidStressArray, ),
    'cpenta_stress_crm' : (RandomSolidStressArray, ),
    'cpenta_stress_psd' : (RandomSolidStressArray, ),
    'cpenta_stress_no' : (RandomSolidStressArray, ),
    'cpenta_stress_rms' : (RandomSolidStressArray, ),

    'chexa_stress' : (RealSolidStressArray, ComplexSolidStressArray),
    'chexa_stress_ato' : (RandomSolidStressArray, ),
    'chexa_stress_crm' : (RandomSolidStressArray, ),
    'chexa_stress_psd' : (RandomSolidStressArray, ),
    'chexa_stress_no' : (RandomSolidStressArray, ),
    'chexa_stress_rms' : (RandomSolidStressArray, ),

    'ctetra_strain' : (RealSolidStrainArray, ComplexSolidStrainArray),
    'ctetra_strain_ato' : (RandomSolidStrainArray, ),
    'ctetra_strain_crm' : (RandomSolidStrainArray, ),
    'ctetra_strain_psd' : (RandomSolidStrainArray, ),
    'ctetra_strain_no' : (RandomSolidStrainArray, ),
    'ctetra_strain_rms' : (RandomSolidStrainArray, ),

    'cpenta_strain' : (RealSolidStrainArray, ComplexSolidStrainArray),
    'cpenta_strain_ato' : (RandomSolidStrainArray, ),
    'cpenta_strain_crm' : (RandomSolidStrainArray, ),
    'cpenta_strain_psd' : (RandomSolidStrainArray, ),
    'cpenta_strain_no' : (RandomSolidStrainArray, ),
    'cpenta_strain_rms' : (RandomSolidStrainArray, ),

    'chexa_strain' : (RealSolidStrainArray, ComplexSolidStrainArray),
    'chexa_strain_ato' : (RandomSolidStrainArray, ),
    'chexa_strain_crm' : (RandomSolidStrainArray, ),
    'chexa_strain_psd' : (RandomSolidStrainArray, ),
    'chexa_strain_no' : (RandomSolidStrainArray, ),
    'chexa_strain_rms' : (RandomSolidStrainArray, ),

    'grid_point_forces' : (RealGridPointForcesArray, ComplexGridPointForcesArray),

    'cquad8_force' : (RealPlateBilinearForceArray, ComplexPlate2ForceArray),
    'cquadr_force' : (RealPlateBilinearForceArray, ComplexPlate2ForceArray),
    'ctria6_force' : (RealPlateBilinearForceArray, ComplexPlate2ForceArray),
    'ctriar_force' : (RealPlateBilinearForceArray, ComplexPlate2ForceArray),

    'ctria3_force' : (RealPlateForceArray, ComplexPlateForceArray),
    'ctria3_force_ato' : (RealPlateForceArray, ),
    'ctria3_force_crm' : (RealPlateForceArray, ),
    'ctria3_force_psd' : (RealPlateForceArray, ),
    'ctria3_force_no' : (RealPlateForceArray, ),
    'ctria3_force_rms' : (RealPlateForceArray, ),

    'cquad4_force' : (RealPlateForceArray, RealPlateBilinearForceArray,
                      ComplexPlateForceArray, ComplexPlate2ForceArray),
    'cquad4_force_ato' : (RealPlateForceArray, ),
    'cquad4_force_crm' : (RealPlateForceArray, ),
    'cquad4_force_psd' : (RealPlateForceArray, ),
    'cquad4_force_no' : (RealPlateForceArray, ),
    'cquad4_force_rms' : (RealPlateForceArray, ),

    'cgap_force' : (RealCGapForceArray, None),

    'cbend_force' : (RealBendForceArray, ComplexCBendForceArray),
    'cconeax_force' : (RealConeAxForceArray, None),

    'cbush_stress' : (RealBushStressArray, ComplexCBushStressArray),
    'cbush_stress_ato' : (RealBushStressArray, ),
    'cbush_stress_crm' : (RealBushStressArray, ),
    'cbush_stress_psd' : (RealBushStressArray, ),
    'cbush_stress_rms' : (RealBushStressArray, ),
    'cbush_stress_no' : (RealBushStressArray, ),

    'cbush_strain' : (RealBushStrainArray, ComplexCBushStrainArray),
    'cbush_strain_ato' : (RealBushStrainArray, ),
    'cbush_strain_crm' : (RealBushStrainArray, ),
    'cbush_strain_psd' : (RealBushStrainArray, ),
    'cbush_strain_rms' : (RealBushStrainArray, ),
    'cbush_strain_no' : (RealBushStrainArray, ),

    'cbush_force' : (RealCBushForceArray, ComplexCBushForceArray), # ComplexCBushForceArray

    'cbush1d_stress_strain' : (RealBush1DStressArray, ComplexCBush1DStressArray),

    'celas1_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'celas2_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'celas3_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'celas4_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'conrod_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'crod_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'ctube_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'cbush_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cbar_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cbeam_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'ctria3_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'ctria6_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'ctriar_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cquad4_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cquad8_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cquadr_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cshear_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cquadx_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'ctetra_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cpenta_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'chexa_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cdum8_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'ctriax_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'ctriax6_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cgap_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'cbend_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'genel_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'dmig_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'nonlinear_celas1_stress' : (RealNonlinearSpringStressArray, ),
    'nonlinear_celas2_stress' : (RealNonlinearSpringStressArray, ),
    'nonlinear_celas3_stress' : (RealNonlinearSpringStressArray, ),
    'nonlinear_celas4_stress' : (RealNonlinearSpringStressArray, ),
    'nonlinear_conrod_stress' : (RealNonlinearRodArray, ),
    'nonlinear_crod_stress' :  (RealNonlinearRodArray, ),
    'nonlinear_ctube_stress' : (RealNonlinearRodArray, ),
    'nonlinear_cgap_stress' : (NonlinearGapStressArray,),
    'nonlinear_cbeam_stress' : (RealNonlinearBeamStressArray, ),

    'nonlinear_cquad4_stress' : (RealNonlinearPlateArray, ),
    'nonlinear_ctria3_stress' : (RealNonlinearPlateArray, ),

    'ctetra_pressure_force' : (RealSolidPressureForceArray, ComplexSolidPressureForceArray,),
    'cpenta_pressure_force' : (RealSolidPressureForceArray, ComplexSolidPressureForceArray,),
    'chexa_pressure_force' : (RealSolidPressureForceArray, ComplexSolidPressureForceArray,),
    'cbeam_force_vu' : (RealCBeamForceVUArray, ComplexCBeamForceVUArray),

    'grid_point_stresses' : (GridPointStressesArray, ),
    'grid_point_volume_stresses' : (GridPointStressesVolumeArray,),

    'crod_thermal_load' :  (Real1DHeatFluxArray,),
    'ctube_thermal_load' :  (Real1DHeatFluxArray,),
    'conrod_thermal_load' : (Real1DHeatFluxArray,),

    'cbar_thermal_load' :  (Real1DHeatFluxArray,),
    'cbeam_thermal_load' : (Real1DHeatFluxArray,),
    'cbend_thermal_load' : (Real1DHeatFluxArray,),

    'cquad4_thermal_load' : (HeatFlux_2D_3DArray,),
    'cquad8_thermal_load' : (HeatFlux_2D_3DArray,),
    'ctria3_thermal_load' : (HeatFlux_2D_3DArray,),
    'ctria6_thermal_load' : (HeatFlux_2D_3DArray,),
    'ctriax6_thermal_load' : (HeatFlux_2D_3DArray,),

    'ctetra_thermal_load' : (HeatFlux_2D_3DArray,),
    'cpenta_thermal_load' : (HeatFlux_2D_3DArray,),
    'chexa_thermal_load' : (HeatFlux_2D_3DArray,),

    'chbdye_thermal_load' :  (RealChbdyHeatFluxArray,),
    'chbdyp_thermal_load' : (RealChbdyHeatFluxArray,),
    'chbdyg_thermal_load' : (RealChbdyHeatFluxArray,),

    'vu_tria_force' : (RealForceVU2DArray,),
    'vu_quad_force' : (RealForceVU2DArray,),

    'temperatures' : (RealTemperatureArray,),
    'thermal_gradient_and_flux' : (RealTemperatureGradientAndFluxArray,),
    'thermal_load_vectors' : (RealTemperatureVectorArray,),
}

TABLE_OBJ_KEYS = list(TABLE_OBJ_MAP.keys())

def _load_eigenvalue(h5_result, log):
    """Loads a RealEigenvalue"""
    class_name = _cast(h5_result.get('class_name'))
    title = ''
    nmodes = _cast(h5_result.get('nmodes'))
    if class_name == 'RealEigenvalues':
        obj = RealEigenvalues(title, nmodes=nmodes)
    elif class_name == 'ComplexEigenvalues':
        obj = ComplexEigenvalues(title, nmodes)
    elif class_name == 'BucklingEigenvalues':
        obj = BucklingEigenvalues(title, nmodes=nmodes)
    else:
        log.warning('  %r is not supported...skipping' % class_name)
        return None

    assert obj.class_name == class_name, 'class_name=%r selected; should be %r' % (obj.class_name, class_name)
    keys_to_skip = ['class_name', 'is_complex', 'is_real']
    for key in h5_result.keys():
        if key in keys_to_skip:
            continue
        else:
            datai = _cast(h5_result.get(key))
            assert not isinstance(datai, binary_type), key
            setattr(obj, key, datai)
    return obj

def _load_table(result_name, h5_result, objs, log, debug=False):# real_obj, complex_obj
    """loads a RealEigenvectorArray/ComplexEigenvectorArray"""
    is_real = _cast(h5_result.get('is_real'))
    #is_complex = _cast(h5_result.get('is_complex'))
    nonlinear_factor = _cast(h5_result.get('nonlinear_factor'))
    #is_stress = _cast(h5_result.get('is_stress'))
    #is_strain = _cast(h5_result.get('is_strain'))

    data_names = [name.decode('utf8') for name in _cast(h5_result.get('data_names')).tolist()]
    str_data_names = [data_name + 's' for data_name in data_names]
    data_code = {
        'nonlinear_factor' : nonlinear_factor,
        'sort_bits' : _cast(h5_result.get('sort_bits')).tolist(),
        'sort_method' : _cast(h5_result.get('sort_method')),
        'is_msc' : _cast(h5_result.get('is_msc')),
        'format_code' : _cast(h5_result.get('format_code')),
        'device_code' : _cast(h5_result.get('device_code')),
        'approach_code' : _cast(h5_result.get('approach_code')),
        'analysis_code' : _cast(h5_result.get('analysis_code')),
        'table_code' : _cast(h5_result.get('table_code')),
        'tCode' : _cast(h5_result.get('tCode')),
        'sort_code' : _cast(h5_result.get('sort_code')),
        'thermal' : _cast(h5_result.get('thermal')),
        'subtitle' : _cast(h5_result.get('subtitle')),
        'acoustic_flag' : _cast(h5_result.get('acoustic_flag')),
        'stress_bits' : _cast(h5_result.get('stress_bits')),
        's_code' : _cast(h5_result.get('s_code')),
        'data_names' : data_names,
        'name' : data_names[0],
        'table_name' : _cast(h5_result.get('table_name')),
    }
    for key, value in list(iteritems(data_code)):
        if value is None:
            if key in ['nonlinear_factor']:
                pass
            elif key in ['acoustic_flag', 'stress_bits', 's_code', 'thermal']:
                del data_code[key]
            else:
                log.warning('%s %s' % (key, value))

    is_sort1 = _cast(h5_result.get('is_sort1'))
    isubcase = _cast(h5_result.get('isubcase'))
    dt = nonlinear_factor

    class_name = _cast(h5_result.get('class_name'))
    obj_class = _get_obj_class(objs, class_name, result_name, is_real, log)
    if obj_class is None:
        log.warning('  unhandled result_name=%r class_name=%r...' % (
            result_name, class_name))
        #raise NotImplementedError('  unhandled result_name=%r class_name=%r...' % (
            #result_name, class_name))
        return None

    obj = obj_class(data_code, is_sort1, isubcase, dt)

    if obj.class_name != class_name:
        msg = 'class_name=%r selected; should be %r' % (obj.class_name, class_name)
        raise RuntimeError(msg)
    _apply_hdf5_attributes_to_object(obj, h5_result, result_name, data_code, str_data_names,
                                     debug=debug)
    return obj

def _apply_hdf5_attributes_to_object(obj, h5_result, result_name, data_code, str_data_names,
                                     debug=False):
    """helper method for ``_load_table``"""
    keys_to_skip = [
        'class_name', 'headers', 'is_real', 'is_complex',
        'is_sort1', 'is_sort2', 'table_name_str',
        'is_curvature', 'is_fiber_distance', 'is_max_shear', 'is_von_mises',
        'is_strain', 'is_stress', 'nnodes_per_element']

    #if result_name == 'eigenvectors':
        #debug = True
    for key in h5_result.keys():
        if key in keys_to_skip:
            continue
        elif result_name == 'grid_point_forces' and key in ['element_name']:
            pass
        elif key in str_data_names:
            if debug:  # pragma: no cover
                print('  *****key=%r' % key)
            datai = _cast(h5_result.get(key))
            setattr(obj, key, datai)
            setattr(obj, '_times', datai)
        elif key not in data_code:
            datai = _cast(h5_result.get(key))
            if debug:  # pragma: no cover
                print('  **key=%r' % key)
                if key not in ['data']:
                    print(datai)

            try:
                setattr(obj, key, datai)
            except AttributeError:
                print('key=%s datai=%r' % (key, datai))
                raise
            if PY3:
                assert not isinstance(datai, binary_type), 'key=%r data=%s' % (key, datai)
    return obj

def _get_obj_class(objs, class_name, result_name, unused_is_real, log):
    #if 1:
    #obj_map = {obj.__class_name : obj for obj in objs if obj is not None}
    #obj_map = {obj.__class__.__name__ : obj for obj in objs if obj is not None}

    # does what the two previous lines should do...
    obj_map = {str(obj).split("'")[1].split('.')[-1] : obj for obj in objs if obj is not None}
    try:
        obj_class = obj_map[class_name]
    except KeyError:
        keysi = list(obj_map.keys())
        print(objs)

        print('obj_map:')
        for key, value in iteritems(obj_map):
            print('  %s : %s' % (key, value))

        # if the obj_map is wrong, you probably have an issue in:
        # - get_oes_prefix_postfix
        # - get_oef_prefix_postfix
        # or other similar function
        log.warning('skipping result_name=%r class_name=%r keys=%s' % (
            result_name, class_name, keysi))
        #raise RuntimeError('result_name=%r class_name=%r keys=%s' % (
            #result_name, class_name, keysi))
        return None
    #else:
        #real_obj, complex_obj = objs
        #if is_real:
            #if real_obj is None:
                #log.warning('    skipping real %r...' % result_name)
                #return None
            #obj_class = real_obj
        #else:
            #if complex_obj is None:
                #log.warning('    skipping complex %r...' % result_name)
                #return None
            #obj_class = complex_obj
    return obj_class

def export_op2_to_hdf5_file(hdf5_filename, op2_model):
    """exports an OP2 object to an HDF5 file"""
    #no_sort2_classes = ['RealEigenvalues', 'ComplexEigenvalues', 'BucklingEigenvalues']

    with h5py.File(hdf5_filename, 'w') as hdf5_file:
        op2_model.log.info('starting export_op2_to_hdf5_file of %r' % hdf5_filename)
        export_op2_to_hdf5(hdf5_file, op2_model)

def export_op2_to_hdf5(hdf5_file, op2_model):
    """exports an OP2 object to an HDF5 file object"""
    info_group = hdf5_file.create_group('info')
    info_group.create_dataset('pyNastran_version', data=pyNastran.__version__)
    info_group.create_dataset('nastran_format', data=op2_model._nastran_format)
    #info_group.create_dataset('is_msc', data=self.is_msc)
    #info_group.create_dataset('is_nx', data=self.is_nx)
    #info_group.create_dataset('nastran_version', data=self.is_nx)
    _export_matrices(hdf5_file, op2_model)
    _export_subcases(hdf5_file, op2_model)

def _export_matrices(hdf5_file, op2_model):
    """exports the matrices to HDF5"""
    if len(op2_model.matrices):
        matrix_group = hdf5_file.create_group('matrices')
        for key, matrix in sorted(iteritems(op2_model.matrices)):
            matrixi_group = matrix_group.create_group(b(key))
            if hasattr(matrix, 'export_to_hdf5'):
                matrix.export_to_hdf5(matrixi_group, op2_model.log)
            else:
                op2_model.log.warning('HDF5: key=%r type=%s cannot be exported' % (key, str(type(matrix))))
                #raise NotImplementedError()
                continue

def _export_subcases(hdf5_file, op2_model):
    """exports the subcases to HDF5"""
    subcase_groups = {}
    result_types = op2_model.get_table_types()
    for result_type in result_types:
        result = getattr(op2_model, result_type)
        #if len(result):
            #print(result)

        for key, obj in iteritems(result):
            #class_name = obj.__class__.__name__
            #print('working on %s' % class_name)
            obj.object_attributes()
            subcase_name = 'Subcase=%s' % str(key)
            if subcase_name in subcase_groups:
                subcase_group = subcase_groups[subcase_name]
            else:
                subcase_group = hdf5_file.create_group(subcase_name)
                subcase_groups[subcase_name] = subcase_group

            #if hasattr(obj, 'element_name'):
                #class_name += ': %s' % obj.element_name

            #result_name = result_type + ':' + class_name
            result_name = result_type
            result_group = subcase_group.create_group(result_name)
            obj.export_to_hdf5(result_group, op2_model.log)

def load_op2_from_hdf5(hdf5_filename, combine=True, log=None):
    """loads an hdf5 file into an OP2 object"""
    assert os.path.exists(hdf5_filename), print_bad_path(hdf5_filename)
    model = OP2(log=None)
    model.op2_filename = hdf5_filename

    log.info('hdf5_op2_filename = %r' % hdf5_filename)
    debug = False
    with h5py.File(hdf5_filename, 'r') as h5_file:
        load_op2_from_hdf5_file(model, h5_file, log, debug=debug)
    model.combine_results(combine=combine)
    return model

def load_op2_from_hdf5_file(model, h5_file, log, debug=False):
    """loads an h5 file object into an OP2 object"""
    for key in h5_file.keys():
        if key.startswith('Subcase'):
            h5_subcase = h5_file.get(key)
            log.debug('subcase:')
            for result_name in h5_subcase.keys():
                if result_name == 'eigenvalues':
                    #log.warning('    skipping %r...' % result_name)
                    h5_result = h5_subcase.get(result_name)
                    obj = _load_eigenvalue(h5_result, log=log)
                    if obj is None:
                        continue
                    model.eigenvalues[obj.title] = obj
                    log.debug('  loaded %r' % result_name)
                elif result_name in TABLE_OBJ_KEYS:
                    if debug:
                        log.debug('  %s:' % result_name)
                    objs = TABLE_OBJ_MAP[result_name]
                    #real_obj, complex_obj = objs
                    h5_result = h5_subcase.get(result_name)
                    if objs is None:
                        log.warning('  skipping %s...' % result_name)
                        continue
                    obj = _load_table(result_name, h5_result, objs, log=log, debug=debug)
                    if obj is None:
                        continue

                    # isubcase, analysis_code, sort_method,
                    #  count, ogs, superelement_adaptivity_index, pval_step
                    opt_count = 0
                    ogs = 0
                    superelement_adaptivity_index = ''
                    pval_step = ''
                    # (1, 2, 1, 0, 0, u'')
                    key = (obj.isubcase, obj.analysis_code, obj.sort_method,
                           opt_count, ogs,
                           superelement_adaptivity_index, pval_step)
                    slot = getattr(model, result_name)
                    slot[key] = obj
                    #log.debug('  loaded %r' % result_name)
                else:
                    log.warning('  unhandled %r...' % result_name)
                    #raise NotImplementedError('  unhandled %r...' % result_name)
            #print(h5_subcase.keys())
        elif key == 'info':
            pass
        elif key == 'matrices':
            _read_h5_matrix(h5_file, model, key, log)
        else:
            log.warning('key = %r' % key)
            #raise NotImplementedError('  unhandled %r...' % key)

def _read_h5_matrix(h5_file, unused_model, key, log):
    """reads an hdf5 matrix"""
    log.debug('matrices:')
    h5_matrix_group = h5_file.get(key)
    for matrix_name in h5_matrix_group.keys():
        h5_matrix = h5_matrix_group.get(matrix_name)
        nkeys = len(h5_matrix.keys())
        if not nkeys:
            log.warning('  %s is empty...skipping' % h5_matrix)
        else:
            log.warning('  skipping %r...' % matrix_name)
            #raise NotImplementedError('matrix=%r' % matrix_name)
            #for attr in h5_matrix.keys():
                #print('    attr=%s' % attr)
