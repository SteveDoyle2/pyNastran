from typing import Dict
from pyNastran.op2.op2_interface.hdf5_interface import (
    RealDisplacementArray, RealBeamStressArray, RealCBeamForceArray, RealGridPointForcesArray)

def check_displacement(displacements: Dict[int, RealDisplacementArray]) -> None:
    unused_case1 = +displacements[1] * -displacements[2] + 2 - 2
    unused_case2 = -displacements[1] * 2 - displacements[2] / 2
    displacements[1] += displacements[2]
    displacements[1] *= displacements[2]
    displacements[1] /= displacements[2]

def check_grid_point_forces(grid_point_forces: Dict[int, RealGridPointForcesArray]) -> None:
    unused_case1 = +grid_point_forces[1] * -grid_point_forces[2] + 2 - 2
    unused_case2 = -grid_point_forces[1] * 2 - grid_point_forces[2] / 2
    grid_point_forces[1] += grid_point_forces[2]
    grid_point_forces[1] *= grid_point_forces[2]
    grid_point_forces[1] /= grid_point_forces[2]
    #grid_point_forces = model.grid_point_forces[1] + model.grid_point_forces[2]

def check_cbeam_stress(cbeam_stress: Dict[int, RealBeamStressArray]) -> None:
    unused_case1 = +cbeam_stress[1] * -cbeam_stress[2] + 2 - 2
    unused_case2 = -cbeam_stress[1] * 2 - cbeam_stress[2] / 2
    cbeam_stress[1] += cbeam_stress[2]
    cbeam_stress[1] *= cbeam_stress[2]
    cbeam_stress[1] /= cbeam_stress[2]

def check_cbeam_force(cbeam_force: Dict[int, RealCBeamForceArray]) -> None:
    unused_case1 = +cbeam_force[1] * -cbeam_force[2] + 2 - 2
    unused_case2 = -cbeam_force[1] * 2 - cbeam_force[2] / 2
    cbeam_force[1] += cbeam_force[2]
    cbeam_force[1] *= cbeam_force[2]
    cbeam_force[1] /= cbeam_force[2]
