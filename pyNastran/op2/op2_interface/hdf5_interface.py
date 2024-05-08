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
from typing import Union, Optional, Any
import numpy as np
import h5py

from cpylog import SimpleLogger
import pyNastran
from pyNastran.op2.op2 import OP2

from pyNastran.op2.result_objects.grid_point_weight import GridPointWeight
from pyNastran.op2.tables.lama_eigenvalues.lama_objects import RealEigenvalues, ComplexEigenvalues, BucklingEigenvalues
from pyNastran.op2.tables.oug.oug_displacements import RealDisplacementArray, ComplexDisplacementArray
from pyNastran.op2.tables.oug.oug_velocities import RealVelocityArray, ComplexVelocityArray
from pyNastran.op2.tables.oug.oug_accelerations import RealAccelerationArray, ComplexAccelerationArray
from pyNastran.op2.tables.oug.oug_eigenvectors import RealEigenvectorArray, ComplexEigenvectorArray
from pyNastran.op2.tables.oug.oug_temperatures import RealTemperatureArray
from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import RealThermalVelocityVectorArray

from pyNastran.op2.result_objects.contact_traction_and_pressure import RealContactTractionAndPressureArray
from pyNastran.op2.tables.contact.slide_objects import RealSlideDistanceArray
from pyNastran.op2.tables.contact.bolt_objects import RealBoltArray



from pyNastran.op2.tables.opg_appliedLoads.opg_load_vector import (
    RealLoadVectorArray, ComplexLoadVectorArray,
    RealTemperatureVectorArray,
)
from pyNastran.op2.tables.oqg_constraintForces.separation_distance import SeparationDistanceArray
from pyNastran.op2.tables.oqg_constraintForces.oqg_contact_forces import RealContactForcesArray
from pyNastran.op2.tables.oqg_constraintForces.oqg_thermal_gradient_and_flux import (
    #RealTemperatureGradientAndFlux,
    RealTemperatureGradientAndFluxArray)

from pyNastran.op2.tables.opr import RealPressureArray
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
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids_nx import RealSolidStrainArrayNx, RealSolidStressArrayNx
from pyNastran.op2.tables.oes_stressStrain.real.oes_solids_composite_nx import RealSolidCompositeStressArray, RealSolidCompositeStrainArray

from pyNastran.op2.tables.oes_stressStrain.real.oes_springs import RealSpringStressArray, RealSpringStrainArray, RealNonlinearSpringStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates import RealCompositePlateStressArray, RealCompositePlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush import RealBushStressArray, RealBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bush1d import RealBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_gap import NonlinearGapStressArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_triax import RealTriaxStressArray #, RealTriaxStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_bend import RealBendStressArray, RealBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_weld import RealWeldStressArray, RealWeldStrainArray
from pyNastran.op2.tables.oes_stressStrain.real.oes_fast import RealFastStressArray, RealFastStrainArray

from pyNastran.op2.tables.oes_stressStrain.real.oes_composite_plates_strength_ratio import RealCompositePlateStressStrengthRatioArray # , RealCompositePlateStrainStrengthRatioArray

from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_rod import RealNonlinearRodArray
from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear import RealNonlinearPlateArray, RealNonlinearSolidArray
from pyNastran.op2.tables.oes_stressStrain.oes_hyperelastic import HyperelasticQuadArray

from pyNastran.op2.tables.oes_stressStrain.complex.oes_bars import ComplexBarStressArray, ComplexBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_beams import ComplexBeamStressArray, ComplexBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush import ComplexCBushStressArray, ComplexCBushStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bush1d import ComplexCBush1DStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates import ComplexPlateStressArray, ComplexPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_composite_plates import ComplexLayeredCompositeStressArray, ComplexLayeredCompositeStrainArray # ComplexLayeredCompositesArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_plates_vm import ComplexPlateVMStressArray, ComplexPlateVMStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_triax import ComplexTriaxStressArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_rods import ComplexRodStressArray, ComplexRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_shear import ComplexShearStressArray, ComplexShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_solids import ComplexSolidStressArray, ComplexSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_springs import ComplexSpringStressArray, ComplexSpringStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_bend import ComplexBendStressArray, ComplexBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_fast import ComplexFastStressArray, ComplexFastStrainArray
from pyNastran.op2.tables.oes_stressStrain.complex.oes_weld import ComplexWeldStressArray, ComplexWeldStrainArray

from pyNastran.op2.tables.oes_stressStrain.random.oes_rods import RandomRodStressArray, RandomRodStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bars import RandomBarStressArray, RandomBarStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_beams import RandomBeamStressArray, RandomBeamStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_bend import RandomBendStressArray, RandomBendStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates import RandomPlateStressArray, RandomPlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_solids import RandomSolidStressArray, RandomSolidStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_shear import RandomShearStressArray, RandomShearStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_composite_plates import RandomCompositePlateStressArray, RandomCompositePlateStrainArray
from pyNastran.op2.tables.oes_stressStrain.random.oes_plates_vm import RandomPlateVMStressArray, RandomPlateVMStrainArray

from pyNastran.op2.tables.oes_stressStrain.oes_nonlinear_bush import RealNonlinearBushArray


from pyNastran.op2.tables.ogs_grid_point_stresses.ogs_surface_stresses import (
    GridPointSurfaceStressesArray, GridPointStressesVolumeDirectArray, GridPointStressesVolumePrincipalArray,
    GridPointStressesSurfaceDiscontinutiesArray,
    # strains
    GridPointSurfaceStrainsArray, GridPointStrainsVolumeDirectArray, GridPointStrainsVolumePrincipalArray,
    GridPointStrainsSurfaceDiscontinutiesArray,
)


from pyNastran.op2.tables.oee_energy.oee_objects import RealStrainEnergyArray, ComplexStrainEnergyArray
from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import RealGridPointForcesArray, ComplexGridPointForcesArray
from pyNastran.op2.tables.oef_forces.oef_force_objects import (
    FailureIndicesArray,
    RealBendForceArray, RealCBar100ForceArray, RealCBarForceArray, RealCBeamForceArray,
    RealCBushForceArray, RealCGapForceArray, RealConeAxForceArray,
    RealCShearForceArray, RealDamperForceArray, RealPlateBilinearForceArray,
    RealPlateForceArray, RealRodForceArray, RealSolidPressureForceArray, RealSpringForceArray,
    RealViscForceArray, RealCBearForceArray,
    RealCWeldForceArray, RealCWeldForceArrayMSC,
    RealCFastForceArrayNX, RealCFastForceArrayMSC,
)
from pyNastran.op2.tables.oef_forces.oef_complex_force_objects import (
    ComplexCBarForceArray, ComplexCBeamForceArray,
    ComplexCBendForceArray, ComplexCBushForceArray, ComplexCShearForceArray,
    ComplexDamperForceArray, ComplexPlate2ForceArray, ComplexPlateForceArray,
    ComplexRodForceArray, ComplexSolidPressureForceArray, ComplexSpringForceArray,
    ComplexViscForceArray, ComplexCBearForceArray,
    ComplexCWeldForceArray, ComplexCWeldForceArrayMSC,
)
from pyNastran.op2.tables.oef_forces.oef_thermal_objects import (
    RealChbdyHeatFluxArray, RealConvHeatFluxArray,
    Real1DHeatFluxArray,
    #RealElementTableArray, RealHeatFluxVUArray,
    RealHeatFlux_2D_3DArray,
)
#from pyNastran.op2.tables.oqg_constraintForces.oqg_thermal_gradient_and_flux import RealTemperatureGradientAndFluxArray
from pyNastran.utils import check_path
from pyNastran.op2.result_objects.matrix import Matrix


def _cast(h5_result_attr):
    """converts the h5py type back into the OP2 type"""
    if h5_result_attr is None:
        return None

    if len(h5_result_attr.shape) == 0:
        out = np.array(h5_result_attr).tolist()
        return out
    out = np.array(h5_result_attr)
    #assert not isinstance(out, str), out
    return out

def _cast_str(h5_result_attr, encoding: str) -> list[str]:
    """converts the h5py type back into the OP2 type"""
    if h5_result_attr is None:
        return None

    if len(h5_result_attr.shape) == 0:
        out = np.array(h5_result_attr).tolist()
        if isinstance(out, bytes):
            out = out.decode(encoding)
        return out
        #raise NotImplementedError(h5_result_attr.dtype)
    else:
        out = np.array(h5_result_attr)
    out2 = [outi.decode(encoding) if isinstance(outi, bytes) else outi
            for outi in out]
    return out2

# the data for these keys must be strings
STRING_KEYS = [
    'result_name', 'superelement_adaptivity_index', 'element_name',
    'label', 'pval_step', 'title']

TABLE_OBJ_MAP = {
    'displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'no.displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'ato.displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'crm.displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'psd.displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'rms.displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'abs.displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'nrl.displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'srss.displacements' : (RealDisplacementArray, ComplexDisplacementArray),
    'acoustic.displacements' : (ComplexDisplacementArray, ),

    'velocities' : (RealVelocityArray, ComplexVelocityArray, RealThermalVelocityVectorArray),
    'no.velocities' : (RealVelocityArray, ComplexVelocityArray),
    'ato.velocities' : (RealVelocityArray, ComplexVelocityArray),
    'crm.velocities' : (RealVelocityArray, ComplexVelocityArray),
    'psd.velocities' : (RealVelocityArray, ComplexVelocityArray),
    'rms.velocities' : (RealVelocityArray, ComplexVelocityArray),
    'abs.velocities' : (RealVelocityArray, ComplexVelocityArray),

    'accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'no.accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'ato.accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'crm.accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'psd.accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'rms.accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'abs.accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'nrl.accelerations' : (RealAccelerationArray, ComplexAccelerationArray),
    'srss.accelerations' : (RealAccelerationArray, ComplexAccelerationArray),

    'solution_set.displacements' : (RealDisplacementArray, ComplexDisplacementArray, ),
    'solution_set.velocities' : (RealVelocityArray, ComplexVelocityArray, ),
    'solution_set.accelerations' : (RealAccelerationArray, ComplexAccelerationArray, ),
    'solution_set.eigenvectors' : (RealEigenvectorArray, ),

    'separation_initial' : (SeparationDistanceArray,),
    'separation_final' : (SeparationDistanceArray,),

    'spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'spc_forces_v' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'no.spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'ato.spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'crm.spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'psd.spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'rms.spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'abs.spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'nrl.spc_forces' : (RealSPCForcesArray, ComplexSPCForcesArray),
    'RAREATC.spc_forces' : (RealSPCForcesArray, ),
    'RARCONS.spc_forces' : (RealSPCForcesArray, ),

    'mpc_forces' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'no.mpc_forces' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'ato.mpc_forces' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'crm.mpc_forces' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'psd.mpc_forces' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'rms.mpc_forces' : (RealMPCForcesArray, ComplexMPCForcesArray),
    'ROQGM1.mpc_forces': (RealMPCForcesArray, ),

    'eigenvectors' : (RealEigenvectorArray, ComplexEigenvectorArray),
    'eigenvectors_structure' : (RealEigenvectorArray, ),
    'eigenvectors_fluid' : (RealEigenvectorArray, ComplexEigenvectorArray, ),
    'RADCONS.eigenvectors' : (RealEigenvectorArray, ),
    'RADEATC.eigenvectors' : (RealEigenvectorArray, ),
    'RADEFFM.eigenvectors' : (RealEigenvectorArray, ),
    'ROUGV1.eigenvectors' : (RealEigenvectorArray, ),

    'force_vectors' : (RealForceVectorArray, ),

    'load_vectors' : (RealLoadVectorArray, ComplexLoadVectorArray),
    'load_vectors_v' : (RealLoadVectorArray, ComplexLoadVectorArray),
    'no.load_vectors' : (RealLoadVectorArray, ),
    'ato.load_vectors' : (RealLoadVectorArray, ),
    'crm.load_vectors' : (RealLoadVectorArray, ),
    'psd.load_vectors' : (RealLoadVectorArray, ),
    'rms.load_vectors' : (RealLoadVectorArray, ),

    'no.pressures' : (RealPressureArray, ),
    'ato.pressures' : (RealPressureArray, ),
    'crm.pressures' : (RealPressureArray, ),
    'psd.pressures' : (RealPressureArray, ),
    'rms.pressures' : (RealPressureArray, ),
    #---------------------------------------------------------------------------
    'stress.celas1_stress' : (RealSpringStressArray, ComplexSpringStressArray),
    'ato.celas1_stress' : (RealSpringStressArray, ),
    'crm.celas1_stress' : (RealSpringStressArray, ),
    'psd.celas1_stress' : (RealSpringStressArray, ),
    'rms.celas1_stress' : (RealSpringStressArray, ),
    'no.celas1_stress' : (RealSpringStressArray, ),
    'modal_contribution.celas1_stress' : (RealSpringStressArray, ComplexSpringStressArray, ),

    'stress.celas2_stress' : (RealSpringStressArray, ComplexSpringStressArray),
    'ato.celas2_stress' : (RealSpringStressArray, ),
    'crm.celas2_stress' : (RealSpringStressArray, ),
    'psd.celas2_stress' : (RealSpringStressArray, ),
    'rms.celas2_stress' : (RealSpringStressArray, ),
    'no.celas2_stress' : (RealSpringStressArray, ),
    'modal_contribution.celas2_stress': (RealSpringStressArray, ComplexSpringStressArray, ),

    'stress.celas3_stress' : (RealSpringStressArray, ComplexSpringStressArray),
    'ato.celas3_stress' : (RealSpringStressArray, ),
    'crm.celas3_stress' : (RealSpringStressArray, ),
    'psd.celas3_stress' : (RealSpringStressArray, ),
    'rms.celas3_stress' : (RealSpringStressArray, ),
    'no.celas3_stress' : (RealSpringStressArray, ),
    'modal_contribution.celas3_stress': (RealSpringStressArray, ComplexSpringStressArray, ),

    'stress.celas4_stress' : (RealSpringStressArray, ComplexSpringStressArray),
    'ato.celas4_stress' : (RealSpringStressArray, ),
    'crm.celas4_stress' : (RealSpringStressArray, ),
    'psd.celas4_stress' : (RealSpringStressArray, ),
    'rms.celas4_stress' : (RealSpringStressArray, ),
    'no.celas4_stress' : (RealSpringStressArray, ),
    'modal_contribution.celas4_stress' : (RealSpringStressArray, ComplexSpringStrainArray, ),  # TODO: I think this is real only...

    'strain.celas1_strain' : (RealSpringStrainArray, ComplexSpringStrainArray),
    'ato.celas1_strain' : (RealSpringStrainArray, ),
    'crm.celas1_strain' : (RealSpringStrainArray, ),
    'psd.celas1_strain' : (RealSpringStrainArray, ),
    'rms.celas1_strain' : (RealSpringStrainArray, ),
    'no.celas1_strain' : (RealSpringStrainArray, ),
    'modal_contribution.celas1_strain' : (RealSpringStrainArray, ComplexSpringStrainArray, ),

    'strain.celas2_strain' : (RealSpringStrainArray, ComplexSpringStrainArray),
    'ato.celas2_strain' : (RealSpringStrainArray, ),
    'crm.celas2_strain' : (RealSpringStrainArray, ),
    'psd.celas2_strain' : (RealSpringStrainArray, ),
    'rms.celas2_strain' : (RealSpringStrainArray, ),
    'no.celas2_strain' : (RealSpringStrainArray, ),
    'modal_contribution.celas2_strain' : (RealSpringStrainArray, ComplexSpringStrainArray, ),

    'strain.celas3_strain' : (RealSpringStrainArray, ComplexSpringStrainArray),
    'ato.celas3_strain' : (RealSpringStrainArray, ),
    'crm.celas3_strain' : (RealSpringStrainArray, ),
    'psd.celas3_strain' : (RealSpringStrainArray, ),
    'rms.celas3_strain' : (RealSpringStrainArray, ),
    'no.celas3_strain' : (RealSpringStrainArray, ),
    'modal_contribution.celas3_strain' : (RealSpringStrainArray, ComplexSpringStrainArray, ),

    'strain.celas4_strain' : (RealSpringStrainArray, ComplexSpringStrainArray),
    'ato.celas4_strain' : (RealSpringStrainArray, ),
    'crm.celas4_strain' : (RealSpringStrainArray, ),
    'psd.celas4_strain' : (RealSpringStrainArray, ),
    'rms.celas4_strain' : (RealSpringStrainArray, ),
    'no.celas4_strain' : (RealSpringStrainArray, ),
    'modal_contribution.celas4_strain' : (RealSpringStrainArray, ComplexSpringStrainArray, ),

    'failure_indices.ctria3_composite_force': (FailureIndicesArray, ),
    'failure_indices.ctria6_composite_force': (FailureIndicesArray, ),
    'failure_indices.ctriar_composite_force': (FailureIndicesArray, ),
    'failure_indices.cquad4_composite_force': (FailureIndicesArray, ),
    'failure_indices.cquad8_composite_force': (FailureIndicesArray, ),
    'failure_indices.cquadr_composite_force': (FailureIndicesArray, ),

    'force.celas1_force' : (RealSpringForceArray, ComplexSpringForceArray),
    'ato.celas1_force' : (RealSpringForceArray, ),
    'crm.celas1_force' : (RealSpringForceArray, ),
    'psd.celas1_force' : (RealSpringForceArray, ),
    'rms.celas1_force' : (RealSpringForceArray, ),
    'no.celas1_force' : (RealSpringForceArray, ),

    'force.celas2_force' : (RealSpringForceArray, ComplexSpringForceArray),
    'ato.celas2_force' : (RealSpringForceArray, ),
    'crm.celas2_force' : (RealSpringForceArray, ),
    'psd.celas2_force' : (RealSpringForceArray, ),
    'rms.celas2_force' : (RealSpringForceArray, ),
    'no.celas2_force' : (RealSpringForceArray, ),

    'force.celas3_force' : (RealSpringForceArray, ComplexSpringForceArray),
    'ato.celas3_force' : (RealSpringForceArray, ),
    'crm.celas3_force' : (RealSpringForceArray, ),
    'psd.celas3_force' : (RealSpringForceArray, ),
    'rms.celas3_force' : (RealSpringForceArray, ),
    'no.celas3_force' : (RealSpringForceArray, ),

    'force.celas4_force' : (RealSpringForceArray, ComplexSpringForceArray),
    'ato.celas4_force' : (RealSpringForceArray, ),
    'crm.celas4_force' : (RealSpringForceArray, ),
    'psd.celas4_force' : (RealSpringForceArray, ),
    'rms.celas4_force' : (RealSpringForceArray, ),
    'no.celas4_force' : (RealSpringForceArray, ),

    'force.cdamp1_force' : (RealDamperForceArray, ComplexDamperForceArray),
    'ato.cdamp1_force' : (RealDamperForceArray, ),
    'crm.cdamp1_force' : (RealDamperForceArray, ),
    'psd.cdamp1_force' : (RealDamperForceArray, ),
    'rms.cdamp1_force' : (RealDamperForceArray, ),
    'no.cdamp1_force' : (RealDamperForceArray, ),

    'force.cdamp2_force' : (RealDamperForceArray, ComplexDamperForceArray),
    'ato.cdamp2_force' : (RealDamperForceArray, ),
    'crm.cdamp2_force' : (RealDamperForceArray, ),
    'psd.cdamp2_force' : (RealDamperForceArray, ),
    'rms.cdamp2_force' : (RealDamperForceArray, ),
    'no.cdamp2_force' : (RealDamperForceArray, ),

    'force.cdamp3_force' : (RealDamperForceArray, ComplexDamperForceArray),
    'ato.cdamp3_force' : (RealDamperForceArray, ),
    'crm.cdamp3_force' : (RealDamperForceArray, ),
    'psd.cdamp3_force' : (RealDamperForceArray, ),
    'rms.cdamp3_force' : (RealDamperForceArray, ),
    'no.cdamp3_force' : (RealDamperForceArray, ),

    'force.cdamp4_force' : (RealDamperForceArray, ComplexDamperForceArray),
    'ato.cdamp4_force' : (RealDamperForceArray, ),
    'crm.cdamp4_force' : (RealDamperForceArray, ),
    'psd.cdamp4_force' : (RealDamperForceArray, ),
    'rms.cdamp4_force' : (RealDamperForceArray, ),
    'no.cdamp4_force' : (RealDamperForceArray, ),

    'force.cvisc_force' : (RealViscForceArray, ComplexViscForceArray),
    'ato.cvisc_force' : (RealViscForceArray, ),
    'crm.cvisc_force' : (RealViscForceArray, ),
    'psd.cvisc_force' : (RealViscForceArray, ),
    'rms.cvisc_force' : (RealViscForceArray, ),
    'no.cvisc_force' : (RealViscForceArray, ),

    'stress.crod_stress' : (RealRodStressArray, ComplexRodStressArray),
    'ato.crod_stress' : (RandomRodStressArray, ),
    'crm.crod_stress' : (RandomRodStressArray, ),
    'psd.crod_stress' : (RandomRodStressArray, ),
    'rms.crod_stress' : (RandomRodStressArray, ),
    'no.crod_stress' : (RandomRodStressArray, ),
    'modal_contribution.crod_stress' : (RealRodStressArray, ComplexRodStressArray, ),

    'stress.conrod_stress' : (RealRodStressArray, ComplexRodStressArray),
    'ato.conrod_stress' : (RandomRodStressArray, ),
    'crm.conrod_stress' : (RandomRodStressArray, ),
    'psd.conrod_stress' : (RandomRodStressArray, ),
    'rms.conrod_stress' : (RandomRodStressArray, ),
    'no.conrod_stress' : (RandomRodStressArray, ),
    'modal_contribution.conrod_stress' : (RealRodStressArray, ComplexRodStressArray, ),

    'stress.ctube_stress' : (RealRodStressArray, ComplexRodStressArray),
    'ato.ctube_stress' : (RandomRodStressArray, ),
    'crm.ctube_stress' : (RandomRodStressArray, ),
    'psd.ctube_stress' : (RandomRodStressArray, ),
    'rms.ctube_stress' : (RandomRodStressArray, ),
    'no.ctube_stress' : (RandomRodStressArray, ),
    'modal_contribution.ctube_stress' : (RealRodStressArray, ComplexRodStressArray, ),

    'strain.crod_strain' : (RealRodStrainArray, ComplexRodStrainArray),
    'ato.crod_strain' : (RandomRodStrainArray, ),
    'crm.crod_strain' : (RandomRodStrainArray, ),
    'psd.crod_strain' : (RandomRodStrainArray, ),
    'rms.crod_strain' : (RandomRodStrainArray, ),
    'no.crod_strain' : (RandomRodStrainArray, ),
    'modal_contribution.crod_strain' : (RealRodStrainArray, ComplexRodStrainArray, ),

    'strain.conrod_strain' : (RealRodStrainArray, ComplexRodStrainArray),
    'ato.conrod_strain' : (RandomRodStrainArray, ),
    'crm.conrod_strain' : (RandomRodStrainArray, ),
    'psd.conrod_strain' : (RandomRodStrainArray, ),
    'rms.conrod_strain' : (RandomRodStrainArray, ),
    'no.conrod_strain' : (RandomRodStrainArray, ),
    'modal_contribution.conrod_strain' : (RealRodStrainArray, ComplexRodStrainArray, ),

    'strain.ctube_strain' : (RealRodStrainArray, ComplexRodStrainArray),
    'ato.ctube_strain' : (RandomRodStrainArray, ),
    'crm.ctube_strain' : (RandomRodStrainArray, ),
    'psd.ctube_strain' : (RandomRodStrainArray, ),
    'rms.ctube_strain' : (RandomRodStrainArray, ),
    'no.ctube_strain' : (RandomRodStrainArray, ),
    'modal_contribution.ctube_strain' : (RealRodStrainArray, ComplexRodStrainArray, ),

    'force.crod_force' : (RealRodForceArray, ComplexRodForceArray),
    'ato.crod_force' : (RealRodForceArray, ),
    'crm.crod_force' : (RealRodForceArray, ),
    'psd.crod_force' : (RealRodForceArray, ),
    'rms.crod_force' : (RealRodForceArray, ),
    'no.crod_force' : (RealRodForceArray, ),

    'force.conrod_force' : (RealRodForceArray, ComplexRodForceArray),
    'ato.conrod_force' : (RealRodForceArray, ),
    'crm.conrod_force' : (RealRodForceArray, ),
    'psd.conrod_force' : (RealRodForceArray, ),
    'rms.conrod_force' : (RealRodForceArray, ),
    'no.conrod_force' : (RealRodForceArray, ),

    'force.ctube_force' : (RealRodForceArray, ComplexRodForceArray),
    'ato.ctube_force' : (RealRodForceArray, ),
    'crm.ctube_force' : (RealRodForceArray, ),
    'psd.ctube_force' : (RealRodForceArray, ),
    'rms.ctube_force' : (RealRodForceArray, ),
    'no.ctube_force' : (RealRodForceArray, ),

    'stress.cbar_stress' : (RealBarStressArray, ComplexBarStressArray),
    'ato.cbar_stress' : (RandomBarStressArray, ),
    'crm.cbar_stress' : (RandomBarStressArray, ),
    'psd.cbar_stress' : (RandomBarStressArray, ),
    'rms.cbar_stress' : (RandomBarStressArray, ),
    'no.cbar_stress' : (RandomBarStressArray, ),
    'modal_contribution.cbar_stress' : (RealBarStressArray, ComplexBarStressArray, ),
    'RASCONS.cbar_stress': (RealBarStressArray, ),
    'RASEATC.cbar_stress': (RealBarStressArray, ),

    'strain.cbar_strain' : (RealBarStrainArray, ComplexBarStrainArray),
    'ato.cbar_strain' : (RandomBarStrainArray, ),
    'crm.cbar_strain' : (RandomBarStrainArray, ),
    'psd.cbar_strain' : (RandomBarStrainArray, ),
    'rms.cbar_strain' : (RandomBarStrainArray, ),
    'no.cbar_strain' : (RandomBarStrainArray, ),
    'modal_contribution.cbar_strain' : (RealBarStrainArray, ComplexBarStrainArray, ),

    'force.cbar_force' : (RealCBarForceArray, RealCBar100ForceArray, ComplexCBarForceArray),
    'abs.cbar_force' : (RealCBarForceArray, ),
    'nrl.cbar_force' : (RealCBarForceArray, ),
    'srss.cbar_force' : (RealCBarForceArray, ),
    'ato.cbar_force' : (RealCBarForceArray, ),
    'crm.cbar_force' : (RealCBarForceArray, ),
    'psd.cbar_force' : (RealCBarForceArray, ),
    'rms.cbar_force' : (RealCBarForceArray, ),
    'no.cbar_force' : (RealCBarForceArray, ),

    'force.cweld_force': (RealCWeldForceArray, RealCWeldForceArrayMSC, ComplexCWeldForceArray, ComplexCWeldForceArrayMSC),
    'force.cfast_force': (RealCFastForceArrayNX, RealCFastForceArrayMSC, ),
    'force.cbear_force': (RealCBearForceArray, ComplexCBearForceArray, ),
    'stress.cfast_stress' : (RealFastStressArray, ComplexFastStressArray, ),
    'strain.cfast_strain' : (RealFastStrainArray, ComplexFastStrainArray, ),

    'stress.cweld_stress' : (RealWeldStressArray, ),

    # lines
    'nrl.cbar_force' : (RealCBarForceArray, ),
    'RAFCONS.cbar_force' : (RealCBarForceArray, ),
    'RAFEATC.cbar_force' : (RealCBarForceArray, ),

    'stress.cbar_stress_10nodes' : (RealBar10NodesStressArray, ),
    'strain.cbar_strain_10nodes' : (RealBar10NodesStrainArray, ),

    'stress.cbeam_stress' : (RealBeamStressArray, ComplexBeamStressArray),
    'ato.cbeam_stress' : (RandomBeamStressArray, ),
    'crm.cbeam_stress' : (RandomBeamStressArray, ),
    'psd.cbeam_stress' : (RandomBeamStressArray, ),
    'rms.cbeam_stress' : (RandomBeamStressArray, ),
    'no.cbeam_stress' : (RandomBeamStressArray, ),
    'modal_contribution.cbeam_stress' : (RealBeamStressArray, ComplexBeamStressArray, ),

    'strain.cbeam_strain' : (RealBeamStrainArray, ComplexBeamStrainArray),
    'ato.cbeam_strain' : (RandomBeamStrainArray, ),
    'crm.cbeam_strain' : (RandomBeamStrainArray, ),
    'psd.cbeam_strain' : (RandomBeamStrainArray, ),
    'rms.cbeam_strain' : (RandomBeamStrainArray, ),
    'no.cbeam_strain' : (RandomBeamStrainArray, ),
    'modal_contribution.cbeam_strain' : (RealBeamStrainArray, ComplexBeamStrainArray, ),

    'stress.cbend_stress' : (RealBendStressArray, ComplexBendStressArray, ),
    'ato.cbend_stress' : (RandomBendStressArray, ),
    'crm.cbend_stress' : (RandomBendStressArray, ),
    'psd.cbend_stress' : (RandomBendStressArray, ),
    'rms.cbend_stress' : (RandomBendStressArray, ),
    'no.cbend_stress' : (RandomBendStressArray, ),

    'strain.cbend_strain' : (RealBendStrainArray, ComplexBendStrainArray, ),
    'ato.cbend_strain' : (RandomBendStrainArray, ),
    'crm.cbend_strain' : (RandomBendStrainArray, ),
    'psd.cbend_strain' : (RandomBendStrainArray, ),
    'rms.cbend_strain' : (RandomBendStrainArray, ),
    'no.cbend_strain' : (RandomBendStrainArray, ),

    'force.cbeam_force' : (RealCBeamForceArray, ComplexCBeamForceArray),
    'ato.cbeam_force' : (RealCBeamForceArray, ),
    'crm.cbeam_force' : (RealCBeamForceArray, ),
    'psd.cbeam_force' : (RealCBeamForceArray, ),
    'rms.cbeam_force' : (RealCBeamForceArray, ),
    'no.cbeam_force' : (RealCBeamForceArray, ),

    # -------------------------------------
    # plates
    'stress.cquad4_stress' : (RealPlateStressArray, ComplexPlateStressArray,
                       # ComplexPlateVMStressArray,
                       ),
    'srss.cquad4_stress' : (RealPlateStressArray,),
    'ato.cquad4_stress' : (RandomPlateStressArray, ),
    'crm.cquad4_stress' : (RandomPlateStressArray, ),
    'psd.cquad4_stress' : (RandomPlateStressArray, ),
    'rms.cquad4_stress' : (RandomPlateStressArray, ),
    'no.cquad4_stress' : (RandomPlateStressArray, RandomPlateVMStressArray),
    'modal_contribution.cquad4_stress' : (RealPlateStressArray, ComplexPlateStressArray, ComplexPlateVMStressArray),
    'RASCONS.cquad4_stress' : (RealPlateStressArray, ),
    'RASEATC.cquad4_stress' : (RealPlateStressArray, ),

    'stress.ctria3_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'ato.ctria3_stress' : (RandomPlateStressArray, ),
    'crm.ctria3_stress' : (RandomPlateStressArray, ),
    'psd.ctria3_stress' : (RandomPlateStressArray, ),
    'rms.ctria3_stress' : (RandomPlateStressArray, ),
    'no.ctria3_stress' : (RandomPlateStressArray, ),
    'modal_contribution.ctria3_stress' : (RealPlateStressArray, ComplexPlateStressArray, ComplexPlateVMStressArray),
    'RASCONS.ctria3_stress' : (RealPlateStressArray, ),
    'RASEATC.ctria3_stress' : (RealPlateStressArray, ),

    'stress.cquad8_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'ato.cquad8_stress' : (RandomPlateStressArray, ),
    'crm.cquad8_stress' : (RandomPlateStressArray, ),
    'psd.cquad8_stress' : (RandomPlateStressArray, ),
    'rms.cquad8_stress' : (RandomPlateStressArray, ),
    'no.cquad8_stress' : (RandomPlateStressArray, ),
    'modal_contribution.cquad8_stress' : (RealPlateStressArray, ComplexPlateStressArray, ComplexPlateVMStressArray),

    'stress.ctria6_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'ato.ctria6_stress' : (RandomPlateStressArray, ),
    'crm.ctria6_stress' : (RandomPlateStressArray, ),
    'psd.ctria6_stress' : (RandomPlateStressArray, ),
    'rms.ctria6_stress' : (RandomPlateStressArray, ),
    'no.ctria6_stress' : (RandomPlateStressArray, ),
    'modal_contribution.ctria6_stress' : (RealPlateStressArray, ComplexPlateStressArray, ComplexPlateVMStressArray),

    'stress.ctriar_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'ato.ctriar_stress' : (RandomPlateStressArray, ),
    'crm.ctriar_stress' : (RandomPlateStressArray, ),
    'psd.ctriar_stress' : (RandomPlateStressArray, ),
    'rms.ctriar_stress' : (RandomPlateStressArray, ),
    'no.ctriar_stress' : (RandomPlateStressArray, ),

    'stress.cquadr_stress' : (RealPlateStressArray, ComplexPlateStressArray),
    'ato.cquadr_stress' : (RandomPlateStressArray, ),
    'crm.cquadr_stress' : (RandomPlateStressArray, ),
    'psd.cquadr_stress' : (RandomPlateStressArray, ),
    'rms.cquadr_stress' : (RandomPlateStressArray, ),
    'no.cquadr_stress' : (RandomPlateStressArray, ),
    'modal_contribution.cquadr_stress' : (RealPlateStressArray, ),

    'strain.cquad4_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ato.cquad4_strain' : (RandomPlateStrainArray, ),
    'crm.cquad4_strain' : (RandomPlateStrainArray, ),
    'psd.cquad4_strain' : (RandomPlateStrainArray, ),
    'rms.cquad4_strain' : (RandomPlateStrainArray, ),
    'no.cquad4_strain' : (RandomPlateStrainArray, ),
    'modal_contribution.cquad4_strain' : (RealPlateStrainArray, ComplexPlateStrainArray, ComplexPlateVMStrainArray),
    'RAECONS.cquad4_strain' : (RealPlateStrainArray, ),
    'RAEEATC.cquad4_strain' : (RealPlateStrainArray, ),

    'strain.ctria3_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ato.ctria3_strain' : (RandomPlateStrainArray, ),
    'crm.ctria3_strain' : (RandomPlateStrainArray, ),
    'psd.ctria3_strain' : (RandomPlateStrainArray, ),
    'rms.ctria3_strain' : (RandomPlateStrainArray, ),
    'no.ctria3_strain' : (RandomPlateStrainArray, ),
    'RAECONS.ctria3_strain' : (RealPlateStrainArray, ),
    'RAEEATC.ctria3_strain' : (RealPlateStrainArray, ),

    'strain.cquad8_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ato.cquad8_strain' : (RandomPlateStrainArray, ),
    'crm.cquad8_strain' : (RandomPlateStrainArray, ),
    'psd.cquad8_strain' : (RandomPlateStrainArray, ),
    'rms.cquad8_strain' : (RandomPlateStrainArray, ),
    'no.cquad8_strain' : (RandomPlateStrainArray, ),
    'modal_contribution.cquad8_strain' : (RealPlateStrainArray, ComplexPlateStrainArray, ComplexPlateVMStrainArray),

    'strain.ctria6_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ato.ctria6_strain' : (RandomPlateStrainArray, ),
    'crm.ctria6_strain' : (RandomPlateStrainArray, ),
    'psd.ctria6_strain' : (RandomPlateStrainArray, ),
    'rms.ctria6_strain' : (RandomPlateStrainArray, ),
    'no.ctria6_strain' : (RandomPlateStrainArray, ),
    'modal_contribution.ctria6_strain' : (RealPlateStrainArray, ComplexPlateStrainArray, ComplexPlateVMStrainArray),

    'strain.ctriar_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ato.ctriar_strain' : (RandomPlateStrainArray, ),
    'crm.ctriar_strain' : (RandomPlateStrainArray, ),
    'psd.ctriar_strain' : (RandomPlateStrainArray, ),
    'rms.ctriar_strain' : (RandomPlateStrainArray, ),
    'no.ctriar_strain' : (RandomPlateStrainArray, ),
    'modal_contribution.ctria3_strain' : (RealPlateStrainArray, ComplexPlateStrainArray, ComplexPlateVMStrainArray),

    'strain.cquadr_strain' : (RealPlateStrainArray, ComplexPlateStrainArray),
    'ato.cquadr_strain' : (RandomPlateStrainArray, ),
    'crm.cquadr_strain' : (RandomPlateStrainArray, ),
    'psd.cquadr_strain' : (RandomPlateStrainArray, ),
    'rms.cquadr_strain' : (RandomPlateStrainArray, ),
    'no.cquadr_strain' : (RandomPlateStrainArray, ),
    'modal_contribution.cquadr_strain' : (RealPlateStrainArray, ),

    'stress.ctriax_stress' : (RealTriaxStressArray, ComplexTriaxStressArray,),

    'stress.cquad4_composite_stress' : (RealCompositePlateStressArray, ),
    'stress.ctria3_composite_stress' : (RealCompositePlateStressArray, ),
    'stress.ctria6_composite_stress' : (RealCompositePlateStressArray, ),
    'stress.ctriar_composite_stress' : (RealCompositePlateStressArray, ),
    'stress.cquad8_composite_stress' : (RealCompositePlateStressArray, ),
    'stress.cquadr_composite_stress' : (RealCompositePlateStressArray, ),

    'strength_ratio.ctria3_composite_stress': (RealCompositePlateStressStrengthRatioArray, ),
    'strength_ratio.ctria6_composite_stress': (RealCompositePlateStressStrengthRatioArray, ),

    'strength_ratio.cquad4_composite_stress' : (RealCompositePlateStressArray,
                                                RealCompositePlateStressStrengthRatioArray),
    'strength_ratio.cquad8_composite_stress': (RealCompositePlateStressStrengthRatioArray, ),

    'strain.cquad4_composite_strain' : (RealCompositePlateStrainArray, ),
    'strain.ctria3_composite_strain' : (RealCompositePlateStrainArray, ),
    'strain.ctria6_composite_strain' : (RealCompositePlateStrainArray, ),
    'strain.ctriar_composite_strain' : (RealCompositePlateStrainArray, ),
    'strain.cquad8_composite_strain' : (RealCompositePlateStrainArray, ),
    'strain.cquadr_composite_strain' : (RealCompositePlateStrainArray, ),

    #'modal_contribution.cquad4_composite_strain' : (RealCompositePlateStrainArray, ),
    'modal_contribution.cquad4_composite_strain' : (ComplexLayeredCompositeStrainArray, ),
    'modal_contribution.cquad4_composite_stress' : (ComplexLayeredCompositeStressArray, ),

    'modal_contribution.cquadr_composite_strain' : (ComplexLayeredCompositeStrainArray, ),
    'modal_contribution.cquadr_composite_stress' : (ComplexLayeredCompositeStressArray, ),

    'modal_contribution.ctria3_composite_strain' : (ComplexLayeredCompositeStrainArray, ),
    'modal_contribution.ctria3_composite_stress' : (ComplexLayeredCompositeStressArray, ),

    'modal_contribution.ctriar_composite_strain' : (ComplexLayeredCompositeStrainArray, ),
    'modal_contribution.ctriar_composite_stress' : (ComplexLayeredCompositeStressArray, ),

    'RAPCONS.cquad4_composite_stress' : (RealCompositePlateStressArray, ),
    'RAPCONS.ctria3_composite_stress' : (RealCompositePlateStressArray, ),
    'RAPCONS.ctria6_composite_stress' : (RealCompositePlateStressArray, ),
    'RAPCONS.cquad8_composite_stress' : (RealCompositePlateStressArray, ),

    'RAPEATC.cquad4_composite_stress' : (RealCompositePlateStressArray, ),
    'RAPEATC.ctria3_composite_stress' : (RealCompositePlateStressArray, ),
    'RAPEATC.ctria6_composite_stress' : (RealCompositePlateStressArray, ),
    'RAPEATC.cquad8_composite_stress' : (RealCompositePlateStressArray, ),

    #'RAPCONS.cquad4_composite_strain' : (RealCompositePlateStrainArray, ),
    #'RAPCONS.ctria3_composite_strain' : (RealCompositePlateStrainArray, ),
    #'RAPCONS.ctria6_composite_strain' : (RealCompositePlateStrainArray, ),
    #'RAPCONS.cquad8_composite_strain' : (RealCompositePlateStrainArray, ),

    'stress.cshear_stress' : (RealShearStressArray, ComplexShearStressArray),
    'ato.cshear_stress' : (RandomShearStressArray, ),
    'crm.cshear_stress' : (RandomShearStressArray, ),
    'psd.cshear_stress' : (RandomShearStressArray, ),
    'rms.cshear_stress' : (RandomShearStressArray, ),
    'no.cshear_stress' : (RandomShearStressArray, ),
    'modal_contribution.cshear_stress' : (ComplexShearStressArray,),

    'strain.cshear_strain' : (RealShearStrainArray, ComplexShearStrainArray),
    'ato.cshear_strain' : (RandomShearStrainArray, ),
    'crm.cshear_strain' : (RandomShearStrainArray, ),
    'psd.cshear_strain' : (RandomShearStrainArray, ),
    'rms.cshear_strain' : (RandomShearStrainArray, ),
    'no.cshear_strain' : (RandomShearStrainArray, ),

    'force.cshear_force' : (RealCShearForceArray, ComplexCShearForceArray),
    'ato.cshear_force' : (RealCShearForceArray, ),
    'crm.cshear_force' : (RealCShearForceArray, ),
    'psd.cshear_force' : (RealCShearForceArray, ),
    'rms.cshear_force' : (RealCShearForceArray, ),
    'no.cshear_force' : (RealCShearForceArray, ),

    'force.coneax_force' : (RealConeAxForceArray, ),

    # -------------------------------------
    # solids
    'stress.ctetra_stress' : (RealSolidStressArray, RealSolidStressArrayNx, ComplexSolidStressArray),
    'ato.ctetra_stress' : (RandomSolidStressArray, ),
    'crm.ctetra_stress' : (RandomSolidStressArray, ),
    'psd.ctetra_stress' : (RandomSolidStressArray, ),
    'rms.ctetra_stress' : (RandomSolidStressArray, ),
    'no.ctetra_stress' : (RandomSolidStressArray, ),
    'RASCONS.ctetra_stress' : (RealSolidStressArray, ),
    'RASEATC.ctetra_stress' : (RealSolidStressArray, ),

    'stress.cpenta_stress' : (RealSolidStressArray, RealSolidStressArrayNx, ComplexSolidStressArray),
    'stress.cpenta_composite_strain': (RealSolidCompositeStrainArray, ),
    'ato.cpenta_stress' : (RandomSolidStressArray, ),
    'crm.cpenta_stress' : (RandomSolidStressArray, ),
    'psd.cpenta_stress' : (RandomSolidStressArray, ),
    'rms.cpenta_stress' : (RandomSolidStressArray, ),
    'no.cpenta_stress' : (RandomSolidStressArray, ),
    'RASCONS.cpenta_stress' : (RealSolidStressArray, ),
    'RASEATC.cpenta_stress' : (RealSolidStressArray, ),

    'stress.chexa_stress' : (RealSolidStressArray, RealSolidStressArrayNx, ComplexSolidStressArray),
    'stress.chexa_composite_stress': (RealSolidCompositeStressArray, ),
    'ato.chexa_stress' : (RandomSolidStressArray, ),
    'crm.chexa_stress' : (RandomSolidStressArray, ),
    'psd.chexa_stress' : (RandomSolidStressArray, ),
    'rms.chexa_stress' : (RandomSolidStressArray, ),
    'no.chexa_stress' : (RandomSolidStressArray, ),
    'RASCONS.chexa_stress' : (RealSolidStressArray, ),
    'RASEATC.chexa_stress' : (RealSolidStressArray, ),

    'strain.cpyram_strain' : (RealSolidStrainArrayNx, ComplexSolidStrainArray, ),
    'stress.cpyram_stress' : (RealSolidStressArrayNx, ComplexSolidStressArray, ),

    'strain.ctetra_strain' : (RealSolidStrainArray, RealSolidStrainArrayNx, ComplexSolidStrainArray),
    'ato.ctetra_strain' : (RandomSolidStrainArray, ),
    'crm.ctetra_strain' : (RandomSolidStrainArray, ),
    'psd.ctetra_strain' : (RandomSolidStrainArray, ),
    'rms.ctetra_strain' : (RandomSolidStrainArray, ),
    'no.ctetra_strain' : (RandomSolidStrainArray, ),
    'RAECONS.ctetra_strain' : (RealSolidStrainArray, ),
    'RAEEATC.ctetra_strain' : (RealSolidStrainArray, ),

    'strain.cpenta_strain' : (RealSolidStrainArray, RealSolidStrainArrayNx, ComplexSolidStrainArray),
    'strain.cpenta_composite_strain': (RealSolidCompositeStrainArray, ),
    'ato.cpenta_strain' : (RandomSolidStrainArray, ),
    'crm.cpenta_strain' : (RandomSolidStrainArray, ),
    'psd.cpenta_strain' : (RandomSolidStrainArray, ),
    'rms.cpenta_strain' : (RandomSolidStrainArray, ),
    'no.cpenta_strain' : (RandomSolidStrainArray, ),
    'RAECONS.cpenta_strain' : (RealSolidStrainArray, ),
    'RAEEATC.cpenta_strain' : (RealSolidStrainArray, ),

    'strain.chexa_strain' : (RealSolidStrainArray, RealSolidStrainArrayNx, ComplexSolidStrainArray),
    'ato.chexa_strain' : (RandomSolidStrainArray, ),
    'crm.chexa_strain' : (RandomSolidStrainArray, ),
    'psd.chexa_strain' : (RandomSolidStrainArray, ),
    'rms.chexa_strain' : (RandomSolidStrainArray, ),
    'no.chexa_strain' : (RandomSolidStrainArray, ),
    'RAECONS.chexa_strain' : (RealSolidStrainArray, ),
    'RAEEATC.chexa_strain' : (RealSolidStrainArray, ),
    'strain.chexa_composite_strain': (RealSolidCompositeStrainArray, ),
    # -------------------------------------

    'grid_point_forces' : (RealGridPointForcesArray, ComplexGridPointForcesArray),
    #'RAGCONS.grid_point_forces' : (RealGridPointForcesArray, ),
    #'RAGEATC.grid_point_forces' : (RealGridPointForcesArray, ),

    'force.cquad8_force' : (RealPlateBilinearForceArray, ComplexPlate2ForceArray),
    'ato.cquad8_force' : (RealPlateBilinearForceArray, ),
    'crm.cquad8_force' : (RealPlateBilinearForceArray, ),
    'psd.cquad8_force' : (RealPlateBilinearForceArray, ),
    'rms.cquad8_force' : (RealPlateBilinearForceArray, ),
    'no.cquad8_force' : (RealPlateBilinearForceArray, ),

    'force.cquadr_force' : (RealPlateForceArray, RealPlateBilinearForceArray, ComplexPlateForceArray, ComplexPlate2ForceArray),
    'ato.cquadr_force' : (RealPlateBilinearForceArray, ),
    'crm.cquadr_force' : (RealPlateBilinearForceArray, ),
    'psd.cquadr_force' : (RealPlateBilinearForceArray, ),
    'rms.cquadr_force' : (RealPlateBilinearForceArray, ),
    'no.cquadr_force' : (RealPlateBilinearForceArray, ),

    'force.ctria6_force' : (RealPlateBilinearForceArray, ComplexPlate2ForceArray),
    'ato.ctria6_force' : (RealPlateBilinearForceArray, ),
    'crm.ctria6_force' : (RealPlateBilinearForceArray, ),
    'psd.ctria6_force' : (RealPlateBilinearForceArray, ),
    'rms.ctria6_force' : (RealPlateBilinearForceArray, ),
    'no.ctria6_force' : (RealPlateBilinearForceArray, ),

    'force.ctriar_force' : (RealPlateForceArray, RealPlateBilinearForceArray, ComplexPlateForceArray, ComplexPlate2ForceArray),
    'ato.ctriar_force' : (RealPlateBilinearForceArray, ),
    'crm.ctriar_force' : (RealPlateBilinearForceArray, ),
    'psd.ctriar_force' : (RealPlateBilinearForceArray, ),
    'rms.ctriar_force' : (RealPlateBilinearForceArray, ),
    'no.ctriar_force' : (RealPlateBilinearForceArray, ),

    'force.ctria3_force' : (RealPlateForceArray, ComplexPlateForceArray),
    'ato.ctria3_force' : (RealPlateForceArray, ),
    'crm.ctria3_force' : (RealPlateForceArray, ),
    'psd.ctria3_force' : (RealPlateForceArray, ),
    'rms.ctria3_force' : (RealPlateForceArray, ),
    'no.ctria3_force' : (RealPlateForceArray, ),
    'RAFCONS.ctria3_force' : (RealPlateForceArray, ), # ?
    'RAFEATC.ctria3_force' : (RealPlateForceArray, ), # ?

    'force.cquad4_force' : (RealPlateForceArray, RealPlateBilinearForceArray,
                            ComplexPlateForceArray, ComplexPlate2ForceArray),
    'ato.cquad4_force' : (RealPlateForceArray, ),
    'crm.cquad4_force' : (RealPlateForceArray, ),
    'psd.cquad4_force' : (RealPlateForceArray, ),
    'rms.cquad4_force' : (RealPlateForceArray, ),
    'no.cquad4_force' : (RealPlateForceArray, RealPlateBilinearForceArray),
    'srss.cquad4_force' : (RealPlateForceArray, ),
    'RAFCONS.cquad4_force' : (RealPlateBilinearForceArray, ),
    'RAFEATC.cquad4_force' : (RealPlateBilinearForceArray, ),

    'force.cgap_force' : (RealCGapForceArray, ),
    'ato.cgap_force' : (RealCGapForceArray, ),
    'crm.cgap_force' : (RealCGapForceArray, ),
    'psd.cgap_force' : (RealCGapForceArray, ),
    'rms.cgap_force' : (RealCGapForceArray, ),
    'no.cgap_force' : (RealCGapForceArray, ),

    'force.cbend_force' : (RealBendForceArray, ComplexCBendForceArray),
    'ato.cbend_force' : (RealBendForceArray, ),
    'crm.cbend_force' : (RealBendForceArray, ),
    'psd.cbend_force' : (RealBendForceArray, ),
    'rms.cbend_force' : (RealBendForceArray, ),
    'no.cbend_force' : (RealBendForceArray, ),

    'force.cconeax_force' : (RealConeAxForceArray, ),

    'stress.cbush_stress' : (RealBushStressArray, ComplexCBushStressArray),
    'ato.cbush_stress' : (RealBushStressArray, ),
    'crm.cbush_stress' : (RealBushStressArray, ),
    'psd.cbush_stress' : (RealBushStressArray, ),
    'rms.cbush_stress' : (RealBushStressArray, ),
    'no.cbush_stress' : (RealBushStressArray, ),
    'modal_contribution.cbush_stress' : (ComplexCBushStressArray, ),
    'nonlinear_cbush_force_stress_strain' : (RealNonlinearBushArray, ),

    'strain.cbush_strain' : (RealBushStrainArray, ComplexCBushStrainArray),
    'ato.cbush_strain' : (RealBushStrainArray, ),
    'crm.cbush_strain' : (RealBushStrainArray, ),
    'psd.cbush_strain' : (RealBushStrainArray, ),
    'rms.cbush_strain' : (RealBushStrainArray, ),
    'no.cbush_strain' : (RealBushStrainArray, ),
    'modal_contribution.cbush_strain' : (ComplexCBushStrainArray, ),

    'force.cbush_force' : (RealCBushForceArray, ComplexCBushForceArray), # ComplexCBushForceArray
    'ato.cbush_force' : (RealCBushForceArray, ),
    'crm.cbush_force' : (RealCBushForceArray, ),
    'psd.cbush_force' : (RealCBushForceArray, ),
    'rms.cbush_force' : (RealCBushForceArray, ),
    'no.cbush_force' : (RealCBushForceArray, ),
    'RAFCONS.cbush_force' : (RealCBushForceArray, ),
    'RAFEATC.cbush_force' : (RealCBushForceArray, ),

    'stress.cbush1d_stress_strain' : (RealBush1DStressArray, ComplexCBush1DStressArray),
    'strain.ctriax_strain' : (RealTriaxStressArray, ),
    'strain.cweld_strain': (RealDisplacementArray, ),

    # -------------------------------------
    # strain energy
    'strain_energy.celas1_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.celas2_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.celas3_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.celas4_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'strain_energy.conrod_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.crod_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.ctube_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'strain_energy.cbush_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cbar_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cbeam_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'strain_energy.ctria3_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.ctria6_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.ctriar_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cquad4_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cquad8_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cquadr_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cshear_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cquadx_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'strain_energy.ctetra_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cpenta_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.chexa_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cpyram_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'strain_energy.cdum8_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.ctriax_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.ctriax6_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cgap_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cbend_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.genel_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.dmig_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.conm2_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.rbe1_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.rbe3_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cweld_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cfast_strain_energy' : (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cseam_strain_energy': (RealStrainEnergyArray, ComplexStrainEnergyArray),
    'strain_energy.cbeam3_strain_energy': (RealStrainEnergyArray, ComplexStrainEnergyArray),

    'RANCONS.cbush_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cbar_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cbeam_strain_energy' : (RealStrainEnergyArray, ),

    'RANCONS.ctria3_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.ctria6_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.ctriar_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cquad4_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cquad8_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cquadr_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cshear_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cquadx_strain_energy' : (RealStrainEnergyArray, ),

    'RANCONS.ctetra_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cpenta_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.chexa_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cdum8_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.ctriax_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.ctriax6_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cgap_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.cbend_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.genel_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.dmig_strain_energy' : (RealStrainEnergyArray, ),
    'RANCONS.conm2_strain_energy' : (RealStrainEnergyArray, ),

    'RANEATC.cbush_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cbar_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cbeam_strain_energy' : (RealStrainEnergyArray, ),

    'RANEATC.ctria3_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.ctria6_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.ctriar_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cquad4_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cquad8_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cquadr_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cshear_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cquadx_strain_energy' : (RealStrainEnergyArray, ),

    'RANEATC.ctetra_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cpenta_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.chexa_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cdum8_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.ctriax_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.ctriax6_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cgap_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.cbend_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.genel_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.dmig_strain_energy' : (RealStrainEnergyArray, ),
    'RANEATC.conm2_strain_energy' : (RealStrainEnergyArray, ),

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
    'nonlinear_cbush1d_stress_strain' : (RealBush1DStressArray, ComplexCBush1DStressArray),
    'nonlinear_ctetra_stress_strain' : (RealNonlinearSolidArray, ),
    'nonlinear_cpenta_stress_strain' : (RealNonlinearSolidArray, ),
    'nonlinear_chexa_stress_strain' : (RealNonlinearSolidArray, ),

    'stress.hyperelastic_cquad4_stress' : (HyperelasticQuadArray, ),
    'strain.hyperelastic_cquad4_strain' : (HyperelasticQuadArray, ),

    'force.ctetra_pressure_force' : (RealSolidPressureForceArray, ComplexSolidPressureForceArray,),
    'force.cpenta_pressure_force' : (RealSolidPressureForceArray, ComplexSolidPressureForceArray,),
    'force.chexa_pressure_force' : (RealSolidPressureForceArray, ComplexSolidPressureForceArray,),
    'force.cpyram_pressure_force' : (RealSolidPressureForceArray, ComplexSolidPressureForceArray,),

    'grid_point_surface_stresses' : (GridPointSurfaceStressesArray, ),
    'grid_point_stresses_volume_direct' : (GridPointStressesVolumeDirectArray, ),
    'grid_point_stresses_volume_principal' : (GridPointStressesVolumePrincipalArray, ),
    'grid_point_stress_discontinuities' : (GridPointStressesSurfaceDiscontinutiesArray, ),

    'grid_point_surface_strains' : (GridPointSurfaceStrainsArray, ),
    'grid_point_strains_volume_direct' : (GridPointStrainsVolumeDirectArray, ),
    'grid_point_strains_volume_principal' : (GridPointStrainsVolumePrincipalArray, ),
    'grid_point_strain_discontinuities' : (GridPointStrainsSurfaceDiscontinutiesArray, ),

    # ----------------------------------------------------------
    # contact
    'contact_slide_distance': (RealSlideDistanceArray, ),
    'contact_tractions_and_pressure': (RealContactTractionAndPressureArray, ),

    'contact_forces' : (RealContactForcesArray, ),
    'glue_forces' : (RealContactForcesArray, ),
    'glue_contact_slide_distance': (RealSlideDistanceArray, ),
    'bolt_results': (RealBoltArray, ),

    # ----------------------------------------------------------
    'thermal_load.conv_thermal_load' : (RealConvHeatFluxArray,),

    'thermal_load.crod_thermal_load' :  (Real1DHeatFluxArray, ),
    'thermal_load.ctube_thermal_load' :  (Real1DHeatFluxArray, ),
    'thermal_load.conrod_thermal_load' : (Real1DHeatFluxArray, ),

    'thermal_load.cbar_thermal_load' :  (Real1DHeatFluxArray, ),
    'thermal_load.cbeam_thermal_load' : (Real1DHeatFluxArray, ),
    'thermal_load.cbend_thermal_load' : (Real1DHeatFluxArray, ),

    'thermal_load.cquad4_thermal_load' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.cquad8_thermal_load' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.ctria3_thermal_load' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.ctria6_thermal_load' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.ctriax6_thermal_load' : (RealHeatFlux_2D_3DArray, ),

    'thermal_load.ctetra_thermal_load' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.cpenta_thermal_load' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.chexa_thermal_load' : (RealHeatFlux_2D_3DArray, ),

    'thermal_load.chbdye_thermal_load' :  (RealChbdyHeatFluxArray, ),
    'thermal_load.chbdyp_thermal_load' : (RealChbdyHeatFluxArray, ),
    'thermal_load.chbdyg_thermal_load' : (RealChbdyHeatFluxArray, ),

    # ----------------------------------------------------------
    'thermal_load.crod_thermal_load_flux' :  (Real1DHeatFluxArray, ),
    'thermal_load.ctube_thermal_load_flux' :  (Real1DHeatFluxArray, ),
    'thermal_load.conrod_thermal_load_flux' : (Real1DHeatFluxArray, ),

    'thermal_load.cbar_thermal_load_flux' :  (Real1DHeatFluxArray, ),
    'thermal_load.cbeam_thermal_load_flux' : (Real1DHeatFluxArray, ),
    'thermal_load.cbend_thermal_load_flux' : (Real1DHeatFluxArray, ),

    'thermal_load.cquad4_thermal_load_flux' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.cquad8_thermal_load_flux' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.ctria3_thermal_load_flux' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.ctria6_thermal_load_flux' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.ctriax6_thermal_load_flux' : (RealHeatFlux_2D_3DArray, ),

    'thermal_load.ctetra_thermal_load_flux' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.cpenta_thermal_load_flux' : (RealHeatFlux_2D_3DArray, ),
    'thermal_load.chexa_thermal_load_flux' : (RealHeatFlux_2D_3DArray, ),

    'thermal_load.chbdye_thermal_load_flux' :  (RealChbdyHeatFluxArray, ),
    'thermal_load.chbdyp_thermal_load_flux' : (RealChbdyHeatFluxArray, ),
    'thermal_load.chbdyg_thermal_load_flux' : (RealChbdyHeatFluxArray, ),
    # ----------------------------------------------------------

    'temperatures' : (RealTemperatureArray, ),
    'thermal_gradient_and_flux' : (RealTemperatureGradientAndFluxArray, ),
    'thermal_load_vectors' : (RealTemperatureVectorArray, ),
}
keya = ''
for key, data in list(TABLE_OBJ_MAP.items()):
    if key.startswith('stress.'):
       keya = key.replace('stress.', 'stressa.')
       TABLE_OBJ_MAP[keya] = data
del keya

TABLE_OBJ_KEYS = list(TABLE_OBJ_MAP.keys())

def _load_grid_point_weight(h5_result):
    """Loads a GridPointWeight"""
    #obj = GridPointWeight()
    #IQ
    #IS
    #MO
    #Q
    #S
    #approach_code
    #cg
    #label
    #mass
    #reference_point
    #subtitle
    #table_code
    #title
    data = {}
    for key in h5_result.keys():
        value = _cast(h5_result.get(key))
        data[key] = value
    datai = []
    for key in ['reference_point', 'MO', 'S', 'mass', 'cg', 'IS', 'IQ', 'Q']:
        datai.append(data[key])
        del data[key]
    obj = GridPointWeight(*datai, *data)
    return obj

def _load_eigenvalue(h5_result, encoding: str,
                     log: SimpleLogger) -> Union[RealEigenvalues, ComplexEigenvalues, BucklingEigenvalues]:
    """Loads a RealEigenvalue"""
    class_name = _cast_str(h5_result.get('class_name'), encoding)
    table_name = '???'
    title = ''
    nmodes = _cast(h5_result.get('nmodes'))
    if nmodes is None:
        log.warning(f'class_name = {class_name} did not define nmodes')
        nmodes = 0
    if class_name == 'RealEigenvalues':
        obj = RealEigenvalues(title, table_name, nmodes=nmodes)
    elif class_name == 'ComplexEigenvalues':
        obj = ComplexEigenvalues(title, table_name, nmodes)
    elif class_name == 'BucklingEigenvalues':
        obj = BucklingEigenvalues(title, table_name, nmodes=nmodes)
    else:
        log.warning('  %r is not supported...skipping' % class_name)
        return None

    assert obj.class_name == class_name, 'class_name=%r selected; should be %r' % (obj.class_name, class_name)
    keys_to_skip = ['class_name', 'is_complex', 'is_real', 'table_name_str']
    for key in h5_result.keys():
        if key in keys_to_skip:
            continue
        else:
            datai = _cast(h5_result.get(key))
            if isinstance(datai, bytes):
                pass
            elif isinstance(datai, str):
                datai = datai.encode('latin1')
            else:
                assert not isinstance(datai, bytes), key
            setattr(obj, key, datai)
    return obj

def _load_table(result_name, h5_result, objs: tuple[Any], encoding: str,
                log: SimpleLogger, debug: bool=False):# real_obj, complex_obj
    """loads a RealEigenvectorArray/ComplexEigenvectorArray"""
    is_real = _cast(h5_result.get('is_real'))
    #is_complex = _cast(h5_result.get('is_complex'))
    nonlinear_factor = _cast(h5_result.get('nonlinear_factor'))
    #is_stress = _cast(h5_result.get('is_stress'))
    #is_strain = _cast(h5_result.get('is_strain'))

    data_names = _cast_str(h5_result.get('data_names'), encoding)
    str_data_names = [data_name + 's' for data_name in data_names]
    data_code = {
        'load_as_h5' : _cast(h5_result.get('load_as_h5')),
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
        'thermal' : _cast_str(h5_result.get('thermal'), encoding),
        'subtitle' : _cast_str(h5_result.get('subtitle'), encoding),
        'acoustic_flag' : _cast(h5_result.get('acoustic_flag')),
        'stress_bits' : _cast(h5_result.get('stress_bits')),
        's_code' : _cast(h5_result.get('s_code')),
        'data_names' : data_names,
        'name' : data_names[0],
        'table_name' : _cast_str(h5_result.get('table_name'), encoding),
    }
    for key, value in list(data_code.items()):
        if isinstance(value, np.ndarray):
            pass
        elif value in (None, np.nan):
            if key in ['nonlinear_factor']:
                pass
            elif key in ['acoustic_flag', 'stress_bits', 's_code', 'thermal']:
                del data_code[key]
            else:
                log.warning(f'{key} {value}')

    is_sort1 = _cast(h5_result.get('is_sort1'))
    isubcase = _cast(h5_result.get('isubcase'))
    dt = nonlinear_factor

    class_name = _cast_str(h5_result.get('class_name'), encoding)

    obj_class = _get_obj_class(objs, class_name, result_name, log)
    if obj_class is None:
        log.warning('  unhandled result_name=%r class_name=%r...' % (
            result_name, class_name))
        raise NotImplementedError('  unhandled result_name=%r class_name=%r...' % (
            result_name, class_name))
        #return None

    obj = obj_class(data_code, is_sort1, isubcase, dt)

    if obj.class_name != class_name:
        msg = 'class_name=%r selected; should be %r' % (obj.class_name, class_name)
        raise RuntimeError(msg)
    _apply_hdf5_attributes_to_object(obj, h5_result, result_name, data_code, str_data_names,
                                     encoding, debug=debug)
    return obj


def _apply_hdf5_attributes_to_object(obj, h5_result, result_name, data_code, str_data_names,
                                     encoding: str, debug: bool=False):
    """helper method for ``_load_table``"""
    keys_to_skip = [
        'class_name', 'headers', 'is_real', 'is_complex',
        'is_sort1', 'is_sort2', 'table_name_str',
        'is_curvature', 'is_fiber_distance', 'is_max_shear', 'is_von_mises',
        'is_strain', 'is_stress', 'nnodes_per_element', 'has_von_mises',
    ]
    filtered_attrs = obj.object_attributes(keys_to_skip=keys_to_skip, filter_properties=True)

    #if result_name == 'eigenvectors':
        #debug = True
    for key in h5_result.keys():
        if key not in filtered_attrs:
            continue
        elif result_name == 'grid_point_forces' and key in ['element_name']:
            pass
        elif key in str_data_names:
            if debug:  # pragma: no cover
                print('  *****key={key!r}')
            datai = _cast_str(h5_result.get(key), encoding)
            setattr(obj, key, datai)
            setattr(obj, '_times', datai)
        elif key not in data_code:
            datai = _cast(h5_result.get(key))
            if debug:  # pragma: no cover
                print('  **key=%r' % key)
                if key not in ['data']:
                    print(datai)

            if key in STRING_KEYS and isinstance(datai, bytes):
                datai = datai.decode('latin1')

            try:
                setattr(obj, key, datai)
            except AttributeError:
                print(f'obj = {obj}')
                print(f'key={key!r} datai={datai!r}')
                raise
            assert not isinstance(datai, bytes), f'key={key!r} data={datai}'
    return obj

def _get_obj_class(objs: tuple[Any],
                   class_name: str,
                   result_name: str,
                   log: SimpleLogger) -> Any:
    """
    Parameters
    ----------
    objs: tuple[Any]
        the tuple of allowable classes -> (RealDisplacementArray, ComplexDisplacementArray)
    class_name: str
        the name of the class (e.g., 'RealDisplacementArray')
    result_name : str
        the type of result
    log : SimpleLogger
        a python logging object

    """
    #if 1:
    #obj_map = {obj.__class_name : obj for obj in objs if obj is not None}
    #obj_map = {obj.__class__.__name__ : obj for obj in objs if obj is not None}

    # a class is a type until it is instantiated...so we have to do some nonsense
    # there's probably a better way to do this
    #
    #  objs = (<class 'pyNastran.op2.tables.oug.oug_displacements.RealDisplacementArray'>,
    #          <class 'pyNastran.op2.tables.oug.oug_displacements.ComplexDisplacementArray'>,)
    #  type(objs[0]) -> 'type'
    #  str(objs[0]) -> 'pyNastran.op2.tables.oug.oug_displacements.RealDisplacementArray'
    #  split the values
    #  'pyNastran', 'op2', 'tables', 'oug', 'oug_displacements', 'RealDisplacementArray'
    #  take the last field
    #  'RealDisplacementArray'
    #  finally map that to a dictionary
    #  obj_map = {
    #      'RealDisplacementArray' : RealDisplacementArray,
    #      'ComplexDisplacementArray' : ComplexDisplacementArray,
    #
    obj_map = {str(obj).split("'")[1].split('.')[-1] : obj
               for obj in objs if obj is not None}
    try:
        obj_class = obj_map[class_name]
    except KeyError:
        keysi = list(obj_map.keys())
        print(objs)

        print('obj_map:')
        for key, value in obj_map.items():
            print(f'  {key} : {value}')

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

def export_op2_to_hdf5_filename(hdf5_filename: str, op2_model: OP2) -> None:
    """exports an OP2 object to an HDF5 file"""
    #no_sort2_classes = ['RealEigenvalues', 'ComplexEigenvalues', 'BucklingEigenvalues']
    try:
        with h5py.File(hdf5_filename, 'w') as hdf5_file:
            op2_model.log.info(f'starting export_op2_to_hdf5_file of {hdf5_filename!r}')
            export_op2_to_hdf5_file(hdf5_file, op2_model)
    except OSError:
        op2_model.log.error(f'failed to export {hdf5_filename!r}')
        raise

def export_op2_to_hdf5_file(hdf5_file, op2_model: OP2) -> None:
    """exports an OP2 object to an HDF5 file object"""
    assert not isinstance(hdf5_file, str), hdf5_file
    create_info_group(hdf5_file, op2_model)
    export_matrices(hdf5_file, op2_model)
    _export_subcases(hdf5_file, op2_model)

def create_info_group(hdf5_file, op2_model: OP2) -> None:
    """creates the info HDF5 group"""
    info_group = hdf5_file.create_group('info')
    info_group.create_dataset('pyNastran_version', data=pyNastran.__version__)
    info_group.create_dataset('nastran_format', data=op2_model._nastran_format)
    #info_group.create_dataset('is_msc', data=self.is_msc)
    #info_group.create_dataset('is_nx', data=self.is_nx)
    #info_group.create_dataset('nastran_version', data=self.is_nx)

def export_matrices(hdf5_file, op2_model: OP2) -> None:
    """exports the matrices to HDF5"""
    if len(op2_model.matrices):
        matrix_group = hdf5_file.create_group('matrices')
        for key, matrix in sorted(op2_model.matrices.items()):
            matrixi_group = matrix_group.create_group(key.encode('latin-1'))
            if hasattr(matrix, 'export_to_hdf5'):
                matrix.export_to_hdf5(matrixi_group, op2_model.log)
            else:
                msg = 'HDF5: key=%r type=%s cannot be exported' % (key, str(type(matrix)))
                op2_model.log.warning(msg)
                raise NotImplementedError(msg)
                #continue

def _export_subcases(hdf5_file, op2_model):
    """exports the subcases to HDF5"""
    subcase_groups = {}
    result_types = op2_model.get_table_types()
    skip_results = ['params', 'gpdt', 'bgpdt', 'eqexin', 'psds', 'cstm', 'trmbu', 'trmbd']
    for result_type in result_types:
        if result_type in skip_results or result_type.startswith('responses.'):
            #op2_model.log.debug('skipping %s' % result_type)
            continue

        result = op2_model.get_result(result_type)
        if result is None:  # gpdt, eqexin
            continue
        #if len(result):
            #print(result)

        for key, obj in result.items():
            #class_name = obj.__class__.__name__
            #print('working on %s' % class_name)
            obj.object_attributes(filter_properties=True)
            subcase_name = 'Subcase=%s' % str(key)
            if '/' in subcase_name:
                name = obj.class_name
                op2_model.log.warning(f"'/' in titles are not supported by HDF5 for {name}; "
                                      "changing to ';'")
                subcase_name = subcase_name.replace('/', ';')
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
    return load_op2_from_hdf5_filename(hdf5_filename, combine=combine, log=log)

def load_op2_from_hdf5_filename(hdf5_filename: str, combine: bool=True,
                                log: Optional[SimpleLogger]=None):
    """loads an hdf5 file into an OP2 object"""
    check_path(hdf5_filename, 'hdf5_filename')
    model = OP2(log=log)
    model.op2_filename = hdf5_filename

    log.info(f'hdf5_op2_filename = {hdf5_filename!r}')
    debug = False
    with h5py.File(hdf5_filename, 'r') as h5_file:
        load_op2_from_hdf5_file(model, h5_file, log, debug=debug)
    model.combine_results(combine=combine)
    return model

def load_op2_from_hdf5_file(model: OP2, h5_file,
                            log: SimpleLogger, debug=False):
    """loads an h5 file object into an OP2 object"""
    encoding = 'latin1'
    for key in h5_file.keys():
        if key.startswith('Subcase'):
            h5_subcase = h5_file.get(key)
            #log.debug('subcase:')
            for result_name in h5_subcase.keys():
                assert isinstance(result_name, str), f'result_name={result_name}; type={type(result_name)}'

                if result_name in ['eigenvalues', 'eigenvalues_fluid']:
                    #log.warning('    skipping %r...' % result_name)
                    h5_result = h5_subcase.get(result_name)
                    obj = _load_eigenvalue(h5_result, encoding, log=log)
                    if obj is None:
                        continue
                    slot = getattr(model, result_name)  # get model.eigenvalues
                    slot[obj.title] = obj
                    log.debug(f'  loaded {result_name!r}')
                elif result_name == 'grid_point_weight':
                    h5_result = h5_subcase.get(result_name)
                    obj = _load_grid_point_weight(h5_result)
                    model.grid_point_weight[obj.superelement_adaptivity_index] = obj

                elif result_name in TABLE_OBJ_KEYS:
                    if debug:
                        log.debug(f'  {result_name}:')
                    objs = TABLE_OBJ_MAP[result_name]
                    #real_obj, complex_obj = objs
                    h5_result = h5_subcase.get(result_name)
                    if objs is None:
                        log.warning(f'  skipping {result_name}...')
                        continue
                    assert isinstance(objs, tuple), f'check that {result_name!r} is tuple in the above dictionary'
                    obj = _load_table(result_name, h5_result, objs,
                                      encoding, log=log, debug=debug)
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
                    slot = model.get_result(result_name)
                    slot[key] = obj
                    #log.debug('  loaded %r' % result_name)
                else:
                    log.warning(f'  unhandled {result_name!r}...')
                    h5_result = h5_subcase.get(result_name)
                    print(h5_result)
                    if 'class_name' in h5_result.keys():
                        print(h5_result['class_name'])
                    raise NotImplementedError(f'  unhandled {result_name}...')
            #print(h5_subcase.keys())
        elif key == 'info':
            pass
        elif key == 'matrices':
            _read_h5_matrix(h5_file, model, key, encoding, log)
        #else:
            #log.warning('key = %r' % key)
            #raise NotImplementedError('  unhandled %r...' % key)

def _read_h5_matrix(h5_file, model, key,
                    encoding: str, log: SimpleLogger):
    """reads an hdf5 matrix"""
    h5_matrix_group = h5_file.get(key)
    matrix_names = []
    matrix_keys = h5_matrix_group.keys()
    for matrix_name in matrix_keys:
        h5_matrix = h5_matrix_group.get(matrix_name)
        nkeys = len(h5_matrix.keys())
        if not nkeys:
            log.warning(f'  {h5_matrix} is empty...skipping')
        else:
            #log.warning('  skipping %r...' % matrix_name)

            #[u'col', u'data', u'form', u'is_matpool', u'name', u'row', u'shape_str']
            name = _cast_str(h5_matrix.get('name'), encoding)
            form = _cast(h5_matrix.get('form'))
            unused_is_matpool = _cast(h5_matrix.get('is_matpool'))
            matrix_obj = Matrix(name, form)

            #matrix = scipy.sparse.coo_matrix(
                #(real_imag, (GCi, GCj)),
                #shape=(mrows, ncols), dtype=dtype)

            #'is_matpool',
            skip_keys = ['name', 'shape_str', 'dtype_str', 'Real', 'Complex', 'GCi', 'GCj', 'GCi_GCj',
                         'is_real', 'is_complex', 'is_dense', 'is_sparse', 'shape', 'tin']
            for keyi in h5_matrix.keys():
                if keyi in skip_keys:
                    continue
                h5_result_attr = h5_matrix.get(keyi)
                value = _cast(h5_result_attr)
                try:
                    setattr(matrix_obj, keyi, value)
                except AttributeError:
                    print(f'    key={key!r} keyi={keyi!r} value={value!r}')
                    raise
            model.matrices[name] = matrix_obj
            matrix_names.append(matrix_name)

    if len(matrix_keys):
        log.debug(f'matrices: {matrix_names}')
