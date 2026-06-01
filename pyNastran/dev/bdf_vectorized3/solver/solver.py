"""
References:
https://www.oasys-software.com/help/gsa/9.0/GSA_Theory.pdf
https://download.strusoft.com/FEM-Design/inst110x/theory.pdf
https://people.duke.edu/~hpgavin/StructuralDynamics/StructuralElements.pdf

Newton Rhapson
--------------
https://www.youtube.com/watch?v=uXvoN4OleeE&ab_channel=Dr.ClaytonPettit
[K][x] = {F0}
[Kt][dx] = {F0} - {F_xi}
{x} += {dx}
{F0}: is the starting force
{Kt}: tangent stiffness

force convergence:
|({F_xi} - {F0})| / |{F0}| < Ftol; Ftol=0.05%

displacement convergence:
|{dx}| / |{dxi}| < utol; utol=

"""

import os
import copy
from datetime import date
from typing import TextIO, Any

import numpy as np
import scipy as sp
import scipy.sparse as sci_sparse
from scipy.sparse import csc_matrix, lil_matrix, dok_matrix, issparse

from cpylog import SimpleLogger

from pyNastran.utils.solver_utils import get_solver
import pyNastran
from pyNastran.nptyping_interface import (
    NDArrayNbool,
    NDArrayNint,
    NDArrayN2int,
    NDArrayNfloat,
    NDArrayNNfloat,
)
from pyNastran.utils.numpy_utils import integer_types  # , float_types
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase

from pyNastran.f06.f06_writer import make_end
from pyNastran.f06.f06_tables.oload_resultant import Resultant

from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_interface.op2_classes import (
    RealDisplacementArray,
    RealSPCForcesArray,
    RealLoadVectorArray,
    ComplexDisplacementArray,
    ComplexVelocityArray,
    ComplexAccelerationArray,
    RealEigenvalues,
    RealEigenvectorArray,
)
from pyNastran.op2.result_objects.grid_point_weight import make_grid_point_weight
# from pyNastran.bdf.mesh_utils.loads import get_ndof

from .recover.utils import get_f06_op2_pch_set, get_mag_phase_from_options

from .recover.freq_force import recover_force_freq
from .recover.modal_force import recover_force_103
from .recover.static_force import recover_force_101

from .recover.static_stress import recover_stress_101
from .recover.static_strain import recover_strain_101
from .recover.strain_energy import recover_strain_energy_101
from .build_stiffness import build_Kgg, DOF_MAP, Kbb_to_Kgg
from .modal_frequency import get_freq_damping
from .partition import partition_matrix, partition_vector, partition_vector2, partition_vector3
from .utils_statics import save_static_table
from .utils_modes import (
    slice_modal_set, get_real_eigenvalue_method,
    apply_phi_normalization, compute_mass_participation,
)
from .utils_freq import get_frequencies, slice_freq_set
from pyNastran.dev.bdf_vectorized3.bdf_interface.breakdowns import NO_MASS  # , NO_VOLUME


def _apply_grav(model, load, scale, Fb, dof_map, ndof, log):
    """Apply GRAV load (gravity) to Fb using grid-level mass approach.

    Following MYSTRAN's approach: F_i = M_grid_i * acceleration.
    The assembled mass matrix is used to compute forces at each grid.
    """
    from .build_stiffness import _COOAccumulator, Kbb_to_Kgg

    for i_load in range(load.n):
        sid = int(load.load_id[i_load])
        N_vec = load.N[i_load]  # direction vector (3,)
        mag = float(load.scale_factor[i_load])  # scale factor (acceleration magnitude)
        cid = int(load.coord_id[i_load])

        # Acceleration vector in basic frame
        norm_N = np.linalg.norm(N_vec)
        if norm_N == 0.0:
            continue
        accel_basic = mag * N_vec / norm_N

        # Transform from CID to basic if needed
        if cid != 0:
            coord = model.coord
            beta = coord.xyz_to_global_transform[cid]
            accel_basic = beta.T @ accel_basic

        accel_basic *= scale

        log.debug(f"  GRAV: accel_basic={accel_basic}")

        # Build lumped mass matrix and apply F = M * a at each grid
        # Use the assembled Mbb to get per-grid mass (6x6 diagonal blocks)
        # More efficient: iterate grids, get diagonal mass, apply accel
        grid = model.grid
        for nid in grid.node_id:
            nid_int = int(nid)
            i1 = dof_map[(nid_int, 1)]
            # Get translational mass from Mbb diagonal (will be assembled later)
            # Instead, apply accel to the force vector; the mass contribution
            # is handled by assembling F = rho*V*a for each element

        # Grid-level approach: we need the assembled mass matrix.
        # Build a temporary Mbb for this calculation.
        # For efficiency, use the same Mbb that will be built for the analysis.
        # The caller (build_Fb) is called AFTER Mbb is built in SOL 101.
        # However, build_Fb is called before Mbb->Mgg transform in the flow.
        # We must build our own mass here.
        from pyNastran.dev.bdf_vectorized3.bdf import Subcase
        temp_subcase = Subcase(id=0)
        Mbb_temp = build_Mbb(model, temp_subcase, dof_map, ndof, fdtype="float64")

        # Apply F = M * a for each grid (translational DOFs only)
        for nid in grid.node_id:
            nid_int = int(nid)
            i1 = dof_map[(nid_int, 1)]
            # Extract the 3x3 translational mass block for this grid
            m_diag = np.array([
                Mbb_temp[i1, i1],
                Mbb_temp[i1 + 1, i1 + 1],
                Mbb_temp[i1 + 2, i1 + 2],
            ])
            Fb[i1: i1 + 3] += m_diag * accel_basic


def _apply_rforce(model, load, scale, Fb, dof_map, ndof, log):
    """Apply RFORCE load (centrifugal/rotational force) to Fb.

    Following MYSTRAN's approach:
    F_i = M_i * [omega x (omega x r_i) + alpha x r_i]
    where r_i = position of grid i relative to the reference grid.
    """
    from pyNastran.dev.bdf_vectorized3.bdf import Subcase

    for i_load in range(load.n):
        nid_ref = int(load.node_id[i_load])
        cid = int(load.coord_id[i_load])
        a_scale = float(load.scale_factor[i_load])
        r1, r2, r3 = load.r123[i_load]
        method = int(load.method[i_load]) if hasattr(load, 'method') else 1
        racc = float(load.racc[i_load]) if hasattr(load, 'racc') else 0.0

        # Angular velocity vector
        omega = a_scale * np.array([r1, r2, r3])

        # Transform from CID to basic if needed
        if cid != 0:
            coord = model.coord
            beta = coord.xyz_to_global_transform[cid]
            omega = beta.T @ omega

        # Angular acceleration (for tangential component)
        alpha = racc * np.array([r1, r2, r3])
        if cid != 0:
            alpha = beta.T @ alpha

        omega *= scale
        alpha *= scale

        log.debug(f"  RFORCE: omega={omega}, alpha={alpha}, nid_ref={nid_ref}")

        # Reference point position
        grid = model.grid
        all_nids = grid.node_id
        xyz_cid0 = grid.xyz_cid0()

        if nid_ref > 0:
            iref = np.searchsorted(all_nids, nid_ref)
            xyz_ref = xyz_cid0[iref]
        else:
            xyz_ref = np.zeros(3)

        # Build temporary Mbb for force computation
        temp_subcase = Subcase(id=0)
        Mbb_temp = build_Mbb(model, temp_subcase, dof_map, ndof, fdtype="float64")

        # Apply F = M * a_centrifugal at each grid
        for nid, xyz_i in zip(all_nids, xyz_cid0):
            nid_int = int(nid)
            i1 = dof_map[(nid_int, 1)]

            # Relative position
            r_i = xyz_i - xyz_ref

            # Centripetal: omega x (omega x r_i)
            omega_cross_r = np.cross(omega, r_i)
            accel_centripetal = np.cross(omega, omega_cross_r)

            # Tangential: alpha x r_i
            accel_tangential = np.cross(alpha, r_i)

            # Total acceleration (negative for inertial force direction in Nastran)
            accel_total = -(accel_centripetal + accel_tangential)

            # Mass at this grid
            m_diag = np.array([
                Mbb_temp[i1, i1],
                Mbb_temp[i1 + 1, i1 + 1],
                Mbb_temp[i1 + 2, i1 + 2],
            ])
            Fb[i1: i1 + 3] += m_diag * accel_total


def _apply_pload1(model, load, scale, Fb, dof_map, log):
    """Apply PLOAD1 (distributed/concentrated loads on CBAR/CBEAM) to Fb."""
    from pyNastran.dev.bdf_vectorized3.solver.beam import (
        beam_pg_distributed,
        beam_pg_point,
        beam_transform,
    )

    for i_load in range(load.n):
        eid = load.element_id[i_load]
        load_type = str(load.load_type[i_load]).strip()
        scale_type = str(load.scale[i_load]).strip()
        x1, x2 = load.x[i_load]
        p1, p2 = load.pressure[i_load]

        # Find the element (CBAR or CBEAM)
        elem = None
        for e in [model.cbar, model.cbeam]:
            if e.n == 0:
                continue
            idx = np.where(e.element_id == eid)[0]
            if len(idx) > 0:
                elem = e.slice_card_by_index(idx)
                break
        if elem is None:
            log.warning(f"  PLOAD1: element {eid} not found")
            continue

        # Get element properties
        xyz1, xyz2 = elem.get_xyz()
        xyz1i = xyz1[0]
        xyz2i = xyz2[0]
        L = np.linalg.norm(xyz2i - xyz1i)
        v, ihat, yhat, zhat, wa, wb = elem.get_axes(xyz1, xyz2)
        e_g_nus = elem.e_g_nu()
        E_val, G_val, nu_val = e_g_nus[0]
        inertia = elem.inertia()
        I1, I2, I12, J_val = inertia[0]
        area = elem.area()[0]
        k_arr = elem.k()
        k1, k2 = k_arr[0]

        # Convert scale type to actual x positions
        if scale_type in ("FR", "FRPR"):
            xa = x1 * L
            xb = x2 * L
        else:
            xa = x1
            xb = x2

        # Determine load direction in element local coords
        # FX/FY/FZ = basic coord; FXE/FYE/FZE = element coord
        is_element = load_type.endswith("E")
        base_type = load_type.rstrip("E")

        qx = qy = qz = 0.0
        qx_end = qy_end = qz_end = 0.0
        is_force = base_type.startswith("F")

        if is_force:
            direction = base_type[1]  # X, Y, or Z
            if is_element:
                if direction == "X":
                    qx, qx_end = p1, p2
                elif direction == "Y":
                    qy, qy_end = p1, p2
                else:
                    qz, qz_end = p1, p2
            else:
                # Basic coord -> element local via rotation
                T = np.vstack([ihat[0], yhat[0], zhat[0]])
                if direction == "X":
                    q_basic = np.array([1.0, 0.0, 0.0])
                elif direction == "Y":
                    q_basic = np.array([0.0, 1.0, 0.0])
                else:
                    q_basic = np.array([0.0, 0.0, 1.0])
                q_local = T @ q_basic
                qx, qy, qz = p1 * q_local
                qx_end, qy_end, qz_end = p2 * q_local
        else:
            # Moment loads (MX, MY, MZ) — not yet supported
            log.warning(f"  PLOAD1: moment load type {load_type} not yet supported")
            continue

        # Point load vs distributed
        is_point = abs(xa - xb) < 1e-14 * L
        if is_point:
            fe = beam_pg_point(
                qx * scale, qy * scale, qz * scale,
                xa, L, E_val, G_val, area, I1, I2, k1, k2,
            )
        else:
            fe = beam_pg_distributed(
                qx * scale, qy * scale, qz * scale,
                L, E_val, G_val, area, I1, I2, k1, k2,
                x_start=xa, x_end=xb,
                qx_end=qx_end * scale, qy_end=qy_end * scale, qz_end=qz_end * scale,
            )

        # Transform to basic and assemble
        Teb = beam_transform(ihat[0], yhat[0], zhat[0])
        PG = Teb.T @ fe

        nid1, nid2 = elem.nodes[0]
        gi1 = dof_map[(nid1, 1)]
        gi2 = dof_map[(nid2, 1)]
        Fb[gi1:gi1 + 6] += PG[:6]
        Fb[gi2:gi2 + 6] += PG[6:]


class Solver:
    """defines the Nastran knockoff class"""

    def __init__(self, model: BDF):
        self.model = model
        model.setup(run_geom_check=True)

        self.superelement_id = 0
        self.op2 = OP2(log=model.log, mode="nx")
        self.log = model.log

        d = date.today()
        # the date stamp used in the F06
        self.op2.date = (d.month, d.day, d.year)

        # the "as solved" a-set displacement (after AUTOSPC is applied)
        self.xa_ = None
        self.Kaa = None
        self.Kgg = None
        self.aset = None
        self.sset = None

        # User-supplied KGG override (e.g. from Nastran DMAP export).
        # When set, the solver skips internal KGG assembly and uses this instead.
        # Must be (ndof, ndof) dense or sparse array in the same DOF ordering
        # as the model's grid nodes (6 DOF per node, node order from grid.node_id).
        self.Kgg_override = None

        base_name = os.path.splitext(model.bdf_filename)[0]
        self._bdf_filename = base_name + ".solver.bdf"
        self.f06_filename = base_name + ".solver.f06"
        self.op2_filename = base_name + ".solver.op2"
        # print(self.f06_filename)
        # print(self.op2_filename)

    def run(self):
        page_num = 1
        model = self.model
        model.write_bdf(self._bdf_filename)
        sol = model.sol
        solmap = {
            31: self.run_sol_31_craig_bampton,
            101: self.run_sol_101_statics,
            103: self.run_sol_103_modes,
            105: self.run_sol_105_buckling,
            108: self.run_sol_108_direct_frequency,
            111: self.run_sol_111_modal_frequency,
        }
        model.cross_reference()
        self._update_card_count()

        # title = ''
        title = f"pyNastran {pyNastran.__version__}"
        for subcase in model.subcases.values():
            if "TITLE" in subcase:
                title = subcase.get_parameter("TITLE")
                break

        today = None
        page_stamp = self.op2.make_stamp(title, today)  # + '\n'
        with open(self.f06_filename, "w") as f06_file:
            f06_file.write(self.op2.make_f06_header())
            self.op2._write_summary(f06_file, card_count=model.card_count)
            f06_file.write("\n")
            if sol in [31, 101, 103, 105, 107, 108, 109, 111, 112]:
                for subcase_id, subcase in sorted(model.subcases.items()):
                    if subcase_id == 0:
                        continue
                    self.log.debug(f"subcase_id={subcase_id}")
                    # isubcase = subcase.id
                    subtitle = f"SUBCASE {subcase_id}"
                    label = ""
                    if "SUBTITLE" in subcase:
                        subtitle = subcase.get_parameter("SUBTITLE")
                    if "LABEL" in subcase:
                        label = subcase.get_parameter("LABEL")

                    runner = solmap[sol]
                    # print(runner)
                    out, page_num, end_options = runner(
                        subcase,
                        f06_file,
                        page_stamp,
                        title=title,
                        subtitle=subtitle,
                        label=label,
                        page_num=page_num,
                        idtype="int32",
                        fdtype="float64",
                    )
                    del out
            else:
                raise NotImplementedError(sol)
            end_flag = True
            f06_file.write(make_end(end_flag, end_options))
        # Write OP2 once after all subcases (only if there's data)
        op2 = self.op2
        has_results = (
            op2.displacements or op2.spc_forces or op2.load_vectors
            or op2.eigenvectors or op2.eigenvalues
        )
        if has_results:
            op2.write_op2(
                self.op2_filename, post=-1, endian=b"<",
                skips=None, nastran_format="nx")

    def _update_card_count(self) -> None:
        for card_type, values in self.model._type_to_id_map.items():
            self.model.card_count[card_type] = len(values)

    def build_Fb(
        self, xg: NDArrayNfloat, sset_b, dof_map: DOF_MAP, ndof: int, subcase: Subcase
    ) -> NDArrayNfloat:
        model = self.model
        log = model.log
        log.info("starting build_Fb")

        Fb = np.zeros(ndof, dtype="float32")

        has_load = "LOAD" in subcase
        has_temp = "TEMPERATURE(LOAD)" in subcase or "TEMPERATURE(BOTH)" in subcase

        if not has_load and not has_temp:
            return Fb

        if has_load:
            load_id, unused_options = subcase["LOAD"]
            reduced_loads = model.get_reduced_static_load()

            loads = reduced_loads[load_id]
            for scale, load in loads:
                if load.type == "SLOAD":
                    for mag, nid in zip(load.mags, load.nodes):
                        i = dof_map[(nid, 0)]  # TODO: wrong...
                        Fb[i] = mag * scale
                elif load.type == "FORCE":
                    # TODO: wrong because it doesn't handle SPOINTs
                    fxyz_myz = load.sum_forces_moments()
                    fxyz = fxyz_myz[:, :3]
                    nids = load.node_id
                    log.debug(f"  FORCE nids={nids} Fxyz={fxyz}")
                    for fxyzi, nid in zip(fxyz, nids):
                        assert len(fxyzi) == 3
                        fi = dof_map[(nid, 1)]
                        Fb[fi : fi + 3] = fxyzi

                elif load.type == "MOMENT":
                    fxyz_myz = load.sum_forces_moments()
                    mxyz = fxyz_myz[:, 3:]
                    nids = load.node_id
                    log.debug(f"  MOMENT nid={nids} Mxyz={mxyz}")
                    for mxyzi, nid in zip(mxyz, nids):
                        fi = dof_map[(nid, 4)]
                        assert len(mxyzi) == 3
                        Fb[fi : fi + 3] = mxyzi
                elif load.type == "SPCD":
                    for nid, components, enforced in zip(
                        load.nodes, load.components, load.enforced
                    ):
                        for component in str(components):
                            dof = int(component)
                            fi = dof_map[(nid, dof)]
                            xg[fi] = enforced
                            Fb[fi] = np.nan
                            sset_b[fi] = True
                elif load.type == "PLOAD1":
                    _apply_pload1(model, load, scale, Fb, dof_map, log)
                elif load.type == "PLOAD4":
                    from .shells import build_pload4_cquad4, build_pload4_ctria3

                    build_pload4_cquad4(model, Fb, dof_map, load_id)
                    build_pload4_ctria3(model, Fb, dof_map, load_id)
                elif load.type == "GRAV":
                    _apply_grav(model, load, scale, Fb, dof_map, ndof, log)
                elif load.type == "RFORCE":
                    _apply_rforce(model, load, scale, Fb, dof_map, ndof, log)
                else:
                    print(load.get_stats())
                    raise NotImplementedError(load)

        # Thermal loads via TEMPERATURE(LOAD) or TEMPERATURE(BOTH)
        if has_temp:
            if "TEMPERATURE(LOAD)" in subcase:
                temp_load_id, _ = subcase["TEMPERATURE(LOAD)"]
            else:
                temp_load_id, _ = subcase["TEMPERATURE(BOTH)"]
            node_temperatures = self._get_node_temperatures(temp_load_id)
            if node_temperatures:
                log.debug(f"  Thermal load: {len(node_temperatures)} nodes with dT")
                from .shells import build_thermal_load_cquad4, build_thermal_load_ctria3
                from .build_stiffness import build_thermal_load_beam

                build_thermal_load_cquad4(model, Fb, dof_map, node_temperatures)
                build_thermal_load_ctria3(model, Fb, dof_map, node_temperatures)
                build_thermal_load_beam(model, Fb, dof_map, node_temperatures)

        log.info("end of build_Fb")
        return Fb

    def _get_node_temperatures(self, temp_load_id: int) -> dict[int, float]:
        """Build a dict of {node_id: temperature} from TEMP/TEMPD cards.

        Parameters
        ----------
        temp_load_id : int
            Load set ID referencing TEMP/TEMPD cards.

        Returns
        -------
        node_temperatures : dict[int, float]
            Temperature at each grid point.
        """
        model = self.model
        node_temperatures: dict[int, float] = {}

        # Default temperature from TEMPD
        default_temp = None
        if model.tempd.n > 0:
            tempd = model.tempd
            idx = np.where(tempd.load_id == temp_load_id)[0]
            if len(idx) > 0:
                default_temp = float(tempd.temperature[idx[0]])

        # Apply default to all grid points
        if default_temp is not None:
            for nid in model.grid.node_id:
                node_temperatures[int(nid)] = default_temp

        # Override with explicit TEMP cards
        if model.temp.n > 0:
            temp = model.temp
            idx = np.where(temp.load_id == temp_load_id)[0]
            for i in idx:
                inode = temp.inode
                i0, i1 = inode[i]
                for j in range(i0, i1):
                    nid = int(temp.node_id[j])
                    t_val = float(temp.temperature[j])
                    node_temperatures[nid] = t_val

        return node_temperatures

    def build_GMN(
        self, subcase: Subcase, dof_map: DOF_MAP, ndof: int,
        fdtype: str = "float64",
    ) -> tuple[np.ndarray | None, np.ndarray | None]:
        """Build the GMN transformation matrix for MPC/rigid element reduction.

        Handles:
          - Explicit MPC cards (from case control MPC = id)
          - RBE2, RBE3, RBAR, RBAR1, RBE1, RROD (always active)

        GMN maps n-set displacements to g-set displacements:
            {u_g} = [GMN] {u_n}

        For independent DOFs, GMN has a 1 on the diagonal.
        For dependent DOFs, GMN has the interpolation coefficients.

        Parameters
        ----------
        subcase : Subcase
            The subcase (checked for MPC card reference).
        dof_map : DOF_MAP
            Maps (nid, dof) to g-set index.
        ndof : int
            Total number of DOFs in the g-set.
        fdtype : str
            Float dtype.

        Returns
        -------
        GMN : (ndof, n_ndof) sparse matrix or None
            Transformation from n-set to g-set. None if no constraints.
        mset : (n_mset,) int array or None
            Indices of dependent (m-set) DOFs in the g-set. None if no constraints.
        """
        model = self.model
        log = model.log

        # Collect all dependent DOFs from MPC cards and rigid elements
        dependents_list = []

        # --- MPC cards ---
        has_mpc = "MPC" in subcase
        mpc_filtered = None
        if has_mpc:
            mpc_id, unused_options = subcase["MPC"]
            mpc = model.mpc
            if mpc.n > 0:
                mpc_filtered = mpc.slice_card_by_id(mpc_id, sort_ids=True)
                if mpc_filtered.n == 0:
                    mpc_filtered = None

        if mpc_filtered is not None:
            for mpc_idi, (idim0, idim1) in zip(mpc_filtered.mpc_id, mpc_filtered.idim):
                nodes = mpc_filtered.node_id[idim0:idim1]
                components = mpc_filtered.components[idim0:idim1]
                nid_dep = int(nodes[0])
                comp_dep = int(components[0])
                idof_dep = dof_map[(nid_dep, comp_dep)]
                dependents_list.append(idof_dep)

        # --- Rigid elements (RBE2, RBE3, RBAR, RBAR1, RBE1, RROD) ---
        has_rigid = _has_rigid_elements(model)
        rigid_gmn_rows = None
        rigid_m_set = None

        if has_rigid:
            from pyNastran.dev.bdf_vectorized3.mesh_utils.gmn_matrix import assemble_gmn
            try:
                rigid_gmn_rows, rigid_m_set, _ = assemble_gmn(
                    model, dof_map=dof_map, ndof=ndof, apply_cd=True)
            except ValueError:
                rigid_gmn_rows = None
                rigid_m_set = None

        if rigid_m_set is not None:
            for (nid, dof), _ in rigid_m_set.items():
                idof_dep = dof_map[(nid, dof)]
                dependents_list.append(idof_dep)

        if not dependents_list:
            return None, None

        mset = np.unique(np.array(dependents_list, dtype="int32"))

        # n-set = g-set minus m-set
        gset = np.arange(ndof, dtype="int32")
        nset = np.setdiff1d(gset, mset)
        n_ndof = len(nset)

        # Map from g-set index to n-set column index
        g_to_n = np.full(ndof, -1, dtype="int32")
        g_to_n[nset] = np.arange(n_ndof, dtype="int32")

        # Build GMN as sparse (ndof x n_ndof)
        # Independent DOFs: identity mapping
        GMN = dok_matrix((ndof, n_ndof), dtype=fdtype)
        for g_idx, n_col in zip(nset, range(n_ndof)):
            GMN[g_idx, n_col] = 1.0

        # --- Fill from MPC cards ---
        if mpc_filtered is not None:
            for mpc_idi, (idim0, idim1) in zip(mpc_filtered.mpc_id, mpc_filtered.idim):
                coefficients = mpc_filtered.coefficients[idim0:idim1]
                components = mpc_filtered.components[idim0:idim1]
                nodes = mpc_filtered.node_id[idim0:idim1]

                nid_dep = int(nodes[0])
                comp_dep = int(components[0])
                coeff_dep = float(coefficients[0])
                idof_dep = dof_map[(nid_dep, comp_dep)]

                for i in range(1, len(nodes)):
                    nid_ind = int(nodes[i])
                    comp_ind = int(components[i])
                    coeff_ind = float(coefficients[i])
                    if coeff_ind == 0.0:
                        continue
                    idof_ind = dof_map[(nid_ind, comp_ind)]
                    n_col = g_to_n[idof_ind]
                    if n_col < 0:
                        log.warning(
                            f"MPC independent DOF ({nid_ind},{comp_ind}) is also "
                            f"dependent in another constraint — skipping")
                        continue
                    GMN[idof_dep, n_col] += -coeff_ind / coeff_dep

        # --- Fill from rigid elements ---
        if rigid_gmn_rows is not None and rigid_m_set is not None:
            # rigid_gmn_rows is (n_m_rigid, ndof) — maps g-set to m-set
            # Each row i corresponds to m-set DOF i, columns are g-set DOFs
            rigid_dense = rigid_gmn_rows.toarray()
            for (nid, dof), m_row in rigid_m_set.items():
                idof_dep = dof_map[(nid, dof)]
                # For each independent g-set DOF contributing to this m-set DOF
                for g_col in range(ndof):
                    val = rigid_dense[m_row, g_col]
                    if abs(val) < 1e-15:
                        continue
                    n_col = g_to_n[g_col]
                    if n_col < 0:
                        continue
                    GMN[idof_dep, n_col] += val

        GMN_csc = GMN.tocsc()
        n_mpc = mpc_filtered.n if mpc_filtered is not None else 0
        n_rigid = len(rigid_m_set) if rigid_m_set is not None else 0
        log.info(f"MPC/rigid reduction: {len(mset)} dependent DOFs eliminated "
                 f"(MPC={n_mpc}, rigid={n_rigid}, g={ndof} -> n={n_ndof})")
        return GMN_csc, mset

    def run_sol_101_statics(
        self,
        subcase: Subcase,
        f06_file: TextIO,
        page_stamp: str,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        page_num: int = 1,
        idtype: str = "int32",
        fdtype: str = "float64",
    ) -> Any:
        """
        Runs a SOL 101

        Analysis (ASET): This set contains all boundary DOFs of the superelement.
                         It is considered fixed by default.
        Fixed Boundary (BSET): This subset of the ASET contains all fixed boundary DOFs.
        Free Boundary  (CSET): This subset of the ASET contains all free boundary DOFs.

        SOL 101 Sets
        ------------
        b = DOFs fixed during component mode analysis or dynamic reduction.
        c = DOFs that are free during component mode synthesis or dynamic reduction.
        lm = Lagrange multiplier DOFs created by the rigid elements
        r = Reference DOFs used to determine free body motion.
        l = b + c + lm
        t = l + r
        q = Generalized DOFs assigned to component modes and residual vectors
        a = t + q
        #-------------------------
        {s} = {sb} + {sg}  # done
        f =
        f = a + o        unconstrained (free) structural DOFs
        n = f + s        all DOFs not constrained by multipoint constraints
        ne = n + e       all DOFs not constrained by multipoint constraints plus extra DOFs
        m = mp + mr      all DOFs eliminated by multipoint constraints
        g = n + m        all DOFs including scalar DOFs

        Per Mystran
        -----------
        n = g - m        all DOFs not constrained by multipoint constraints
        f = n - s        unconstrained (free structural DOFs)
        a = f - o        ananlysis??? set
        L = a - r

        My Simplification
        -----------------
        g = n + (mp + mr)
        g = f + s + (mp + mr)
        g = f + (sb + sg) + (mp + mr)
        g = (a + o) + (sb + sg) + (mp + mr)
        a + o = g - (sb + sg) + (mp + mr)
        a = g - o - (sb + sg) + (mp + mr)

        l = b + c   # left over dofs
        t = l + r   # total set of DOFs for super
        a = t + q
        d = a + e   # dynamic analysis
        f = a + o
        n = f + s
        g = n + m
        p = g + e
        ps = p + sa
        pa = ps + k
        fe = f + e
        ne = n + e
        fr = f - q - r
        v = o + c + r

        New eqs
        f = n - s
        a = f - o
        t = a - q
        l = t - r
        b = l - c ???
        c = l - b ???
        """
        self.log.debug(f"run_sol_101")
        end_options = [
            "SEMG",  # STIFFNESS AND MASS MATRIX GENERATION STEP
            "SEMR",  # MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES)
            "SEKR",  # STIFFNESS MATRIX REDUCTION STEP
            "SELG",  # LOAD MATRIX GENERATION STEP
            "SELR",  # LOAD MATRIX REDUCTION STEP
        ]
        # basic
        itime = 0
        ntimes = 1  # static
        isubcase = subcase.id
        # -----------------------------------------------------------------------
        model = self.model
        model.setup(run_geom_check=True)

        log = model.log

        dof_map, ps = _get_dof_map(model)

        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        ngrid, ndof_per_grid, ndof = get_ndof(model, subcase)

        gset_b = ps_to_sg_set(ndof, ps)
        if self.Kgg_override is not None:
            Kgg_in = self.Kgg_override
            if issparse(Kgg_in):
                Kgg = Kgg_in[:ndof, :ndof].tocsc()
            else:
                Kgg = csc_matrix(Kgg_in[:ndof, :ndof])
        else:
            Kgg = build_Kgg(model, dof_map, ndof, ngrid, ndof_per_grid, idtype="int32", fdtype=fdtype)
        Mbb = build_Mbb(model, subcase, dof_map, ndof, fdtype=fdtype)
        if Mbb is not None:
            reference_point, MO = grid_point_weight(model, Mbb, dof_map, ndof)
            weight = make_grid_point_weight(
                reference_point,
                MO,
                approach_code=1,
                table_code=13,
                title=title,
                subtitle=subtitle,
                label=label,
                superelement_adaptivity_index="",
            )
            self.op2.grid_point_weight[label] = weight
            page_num = weight.write_f06(f06_file, page_stamp, page_num)

        # ----------------------MPC reduction (g -> n)----------------------
        GMN, mset = self.build_GMN(subcase, dof_map, ndof, fdtype=fdtype)
        self.GMN = GMN
        self.mset = mset

        Mgg = Kbb_to_Kgg(model, Mbb, ngrid, ndof_per_grid, inplace=False) if Mbb is not None else None

        if GMN is not None:
            # Transform to n-set: Knn = GMN^T @ Kgg @ GMN, Mnn = GMN^T @ Mgg @ GMN
            Knn = GMN.T @ Kgg @ GMN
            Mnn = GMN.T @ Mgg @ GMN if Mgg is not None else None

            # AUTOSPC on N-set: detect singular DOFs after MPC reduction
            n_ndof_temp = Knn.shape[0]
            autospc_dofs_n = autospc_n_set(Knn, n_ndof_temp, log)
            if len(autospc_dofs_n) > 0:
                # Add these to the SPC set in the n-set coordinate system
                # They will be picked up during the SPC partition below
                self._autospc_n_dofs = autospc_dofs_n
            else:
                self._autospc_n_dofs = np.array([], dtype="int32")
        else:
            Knn = Kgg
            Mnn = Mgg
            self._autospc_n_dofs = np.array([], dtype="int32")

        gset = np.arange(ndof, dtype=idtype)
        sset, sset_b, xg = _build_xg(model, dof_map, ndof, subcase)
        Fb = self.build_Fb(xg, sset_b, dof_map, ndof, subcase)
        Fg = Fb

        # Save g-set force for oload output before MPC transform
        Fg_gset = np.where(np.isnan(Fg), 0.0, Fg.copy())

        if GMN is not None:
            # Transform force to n-set
            Fn = GMN.T @ Fg
            # Transform enforced displacements to n-set
            # n-set indices: remove m-set DOFs from the indexing
            nset_indices = np.setdiff1d(gset, mset)
            n_ndof = len(nset_indices)

            # Rebuild SPC set and xg in n-set coordinates
            # Map g-set SPC indices to n-set indices
            g_to_n = np.full(ndof, -1, dtype="int32")
            g_to_n[nset_indices] = np.arange(n_ndof, dtype="int32")

            sset_n = []
            for s_idx in sset:
                n_idx = g_to_n[s_idx]
                if n_idx >= 0:
                    sset_n.append(n_idx)
            sset_n = np.array(sset_n, dtype="int32")

            xn = np.full(n_ndof, np.nan, dtype=fdtype)
            for s_idx in sset:
                n_idx = g_to_n[s_idx]
                if n_idx >= 0:
                    xn[n_idx] = xg[s_idx]

            # Also check if any dependent DOF has an SPC — that's a conflict
            for m_idx in mset:
                if sset_b[m_idx]:
                    log.warning(
                        f"DOF {m_idx} is both MPC-dependent and SPC-constrained; "
                        f"SPC on dependent DOF is ignored")

            # Add AUTOSPC DOFs to the SPC set in n-set
            if len(self._autospc_n_dofs) > 0:
                sset_n = np.union1d(sset_n, self._autospc_n_dofs).astype("int32")
                for a_idx in self._autospc_n_dofs:
                    if np.isnan(xn[a_idx]):
                        xn[a_idx] = 0.0

            # Work in n-set from here
            nset_all = np.arange(n_ndof, dtype=idtype)
            sset_b_n = np.zeros(n_ndof, dtype="bool")
            sset_b_n[sset_n] = True

            asetmap = get_aset(model)
            if asetmap:
                aset_g = apply_dof_map_to_set(asetmap, dof_map, idtype=idtype)
                aset = np.array([g_to_n[i] for i in aset_g if g_to_n[i] >= 0],
                                dtype=idtype)
            else:
                aset = np.setdiff1d(nset_all, sset_n)

            sset = sset_n
            xg_orig = xg
            xg = xn
            Fg = Fn
            Kgg_orig = Kgg
            Kgg = Knn
            ndof_solve = n_ndof
        else:
            ndof_solve = ndof
            Kgg_orig = Kgg
            xg_orig = None

            # ----------------------SPC reduction (g -> a)----------------------
            asetmap = get_aset(model)
            if asetmap:
                aset = apply_dof_map_to_set(asetmap, dof_map, idtype=idtype)
            else:
                aset = np.setdiff1d(gset, sset)

        naset = aset.sum() if aset.dtype == bool else len(aset)
        nsset = sset.sum() if sset.dtype == bool else len(sset)
        if naset == 0 and nsset == 0:
            raise RuntimeError("no residual structure found")
        # The a-set and o-set are created in the following ways:
        # 1. If only OMITi entries are present, then the o-set consists
        #    of DOFs listed explicitly on OMITi entries. The remaining
        #    f-set DOFs are placed in the b-set, which is a subset of
        #    the a-set.
        # 2. If ASETi or QSETi entries are present, then the a-set consists
        #    of all DOFs listed on ASETi entries and any entries listing its
        #    subsets, such as QSETi, SUPORTi, CSETi, and BSETi entries.  Any
        #    OMITi entries are redundant. The remaining f-set DOFs are placed
        #    in the o-set.
        # 3. If there are no ASETi, QSETi, or OMITi entries present but there
        #    are SUPORTi, BSETi, or CSETi entries present, then the entire
        #    f-set is placed in the a-set and the o-set is not created.

        Fg_oload = Fg_gset
        page_num = write_oload(
            Fg_oload, dof_map, isubcase, ngrid, ndof_per_grid, f06_file, page_stamp, page_num, log
        )

        # aset - analysis set
        # sset - SPC set
        # print('aset = ', aset)
        # print('sset = ', sset)

        # u1 = Kaa^-1 * (F1k - Kas*u2k)
        abs_xg = np.abs(xg)
        finite_xg = np.any(np.isfinite(abs_xg))
        if finite_xg and np.nanmax(abs_xg) > 0.0:
            self.log.info(f"SPCD found")
            self.log.info(f"  xg = {xg}")
            self.log.info(f"  Fg = {Fg}")
            set0 = xg == 0.0
            set0_ = np.where(set0)

            self.log.info(f"  aset = {aset}")
            if np.any(set0_):
                self.log.info(f"  removing set0_ from sset")
                self.log.info(f"    set0_ = {set0_}")
                self.log.info(f"    sset = {sset}")
                sset = np.setdiff1d(sset, set0_)
                self.log.info(f"    -> sset = {sset}")
            else:
                self.log.info(f"  sset = {sset}")

            # self.log.info(f'  Fa_solve = {Fa_solve}')
            # Fa = Fa.reshape(len(Fa), 1)
            # xa = Fa.reshape(len(xa), 1)
            # print(f'Fa.shape = {Fa.shape}')
            # print(f'xa.shape = {xa.shape}')
            # print(f'Kas.shape = {Kas.shape}')
            # Kas_xa = Ksa @ xa
            # print(f'Kas @ xa.shape = {Kas_xa.shape}')
            # Fa_solve = Fa - Kas @ xa
            # Kspc = Kas[ixa, :][:, ixa]
            # print(Kspc.shape)
        else:
            set0 = sset
            sset = []
            # self.sset = sset
        self.sset = sset
        self.set0 = set0
        self.aset = aset

        Fa, Fs = partition_vector2(Fg, [["a", aset], ["s", sset]])
        del Fs

        xa, xs, x0 = partition_vector3(xg, [["a", aset], ["s", sset], ["0", set0]])
        # self.log.info(f'xg = {xg}')
        del xg
        # self.log.info(f'xa = {xa}')
        # self.log.info(f'xs = {xs}')
        # self.log.info(f'x0 = {x0}')
        del x0

        self.Kgg = Kgg_orig
        K = partition_matrix(Kgg, [("a", aset), ("s", sset), ("0", set0)])
        Kaa = K["aa"]
        Kss = K["ss"]
        # Kas = K['as']
        Ksa = K["sa"]
        # Ks = partition_matrix(Kggs, [['a', aset], ['s', sset], ['0', set0]])
        # Kaas = Ks['aa']
        # assert Kaa.shape == Kaas.shape

        self.Kaa = Kaa

        # M = partition_matrix(Mgg, [['a', aset], ['s', sset]])
        # Maa = M['aa']
        # Kss = K['ss']
        Kas = K["as"]
        Ksa = K["sa"]
        K0a = K["0a"]
        K0s = K["0s"]
        # self.Maa = Maa
        # [Kaa]{xa} + [Kas]{xs} = {Fa}
        # [Ksa]{xa} + [Kss]{xs} = {Fs}

        # {xa} = [Kaa]^-1 * ({Fa} - [Kas]{xs})
        # {Fs} = [Ksa]{xa} + [Kss]{xs}

        # print(Kaa)
        # print(Kas)
        # print(Kss)
        Fa_solve: np.ndarray = Fa
        is_sset = len(sset)
        is_set0 = len(set0)
        is_aset = len(aset)
        if is_sset:
            Fa_solve = Fa - Kas @ xs
            self.log.info(f"  Fa_solve = {Fa_solve}")

        # --- SUPORT / inertia relief for statics ---
        rset_b = get_rset_bool(model, dof_map, ndof)
        has_suport = np.any(rset_b)
        inrel = -1
        if hasattr(model, 'params') and 'INREL' in model.params:
            inrel = model.params['INREL'].values[0]

        if has_suport and inrel == -2 and is_aset:
            # Inertia relief: partition a = l + r, apply inertia relief
            from pyNastran.dev.bdf_vectorized3.mesh_utils.inertia_relief import (
                build_rigid_body_modes, compute_inertia_relief,
            )
            a_indices = np.where(aset)[0]
            r_in_a = rset_b[a_indices]
            lset_local = np.where(~r_in_a)[0]
            rset_local = np.where(r_in_a)[0]
            nr = len(rset_local)
            nl = len(lset_local)
            log.info(
                f"Inertia relief (PARAM,INREL,-2): "
                f"r-set={nr} DOFs, l-set={nl} DOFs")

            # Build D matrix for a-set grid coordinates
            # Map a-set DOF indices back to grid positions
            node_xyz = model.grid.xyz_cid0()
            grid_ids = model.grid.node_id
            D_full = build_rigid_body_modes(grid_ids, node_xyz)

            # Extract D for a-set DOFs only
            D_a = D_full[a_indices, :]

            # Build a-set mass matrix (dense)
            M_a = partition_matrix(
                Mgg, [("a", aset), ("s", sset), ("0", set0)])
            Maa_dense = M_a["aa"]
            if hasattr(Maa_dense, 'toarray'):
                Maa_dense = Maa_dense.toarray()

            # Apply inertia relief
            F_net, a_rigid, F_inertia = compute_inertia_relief(
                Maa_dense, D_a, Fa_solve)
            log.info(f"  rigid body acceleration: {a_rigid}")

            # Partition to l-set and solve
            Kaa_dense = Kaa
            if hasattr(Kaa_dense, 'toarray'):
                Kaa_dense = Kaa_dense.toarray()
            Kll = Kaa_dense[np.ix_(lset_local, lset_local)]
            Fl = F_net[lset_local]

            # Solve Kll @ ul = Fl
            Kll_sp = csc_matrix(Kll)
            xl_, ipositive_l, inegative_l = solve(
                Kll_sp, Fl, np.ones(nl, dtype='bool'), log,
                idtype=idtype)

            # Expand l-set solution back to a-set
            xa[:] = 0.0
            ul = np.zeros(nl, dtype=fdtype)
            ul[ipositive_l] = xl_
            xa[lset_local] = ul
            # r-set displacements are zero (rigid body reference)

            self.xa_ = xa[lset_local][ipositive_l]
            self.Fa_ = Fl[ipositive_l]
            self.inertia_relief = {
                'a_rigid': a_rigid,
                'F_inertia': F_inertia,
                'F_net': F_net,
                'D_a': D_a,
            }
        elif is_aset:
            xa_, ipositive, inegative = solve(Kaa, Fa_solve, aset, log, idtype=idtype)
            Fa_ = Fa[ipositive]

            log.info(f"aset_ = {ipositive}")
            log.info(f"xa_ = {xa_}")
            log.info(f"Fa_ = {Fa_}")

            xa[ipositive] = xa_
            xa[inegative] = 0.0
            self.xa_ = xa_
            self.Fa_ = Fa_
        else:
            self.log.warning("A-set is empty; all DOFs are constrained")
            self.xa_ = []
            self.Fa_ = []

        # Assemble solution in the working DOF set (n-set if GMN, else g-set)
        xn_full = np.full(ndof_solve, np.nan, dtype=fdtype)
        xn_full[aset] = xa
        xn_full[sset] = xs
        xn_full[set0] = 0.0
        Fg[aset] = Fa

        fspc = np.full(ndof_solve, 0.0, dtype=fdtype)

        if is_sset:
            log.info(f"fspc_s recovery")
            log.debug(f"  aset = {aset}")
            log.debug(f"  sset = {sset}")
            log.debug(f"  xa = {xa}")
            log.debug(f"  xs = {xs}")
            fspc_s = Ksa @ xa + Kss @ xs
            log.debug(f"  fspc_s = {fspc_s}")
            log.debug(f"  Ksa @ xa = {Ksa @ xa}")
            log.debug(f"  Kss @ xs = {Kss @ xs}")
            log.debug(f"  sum = {Ksa @ xa + Kss @ xs}")
            Fg[sset] = fspc_s
            fspc[sset] = fspc_s
        if is_set0:
            log.info(f"fspc_0 recovery")
            fspc_0 = K0a @ xa + K0s @ xs
            log.debug(f"  fspc_0 = {fspc_0}")
            Fg[set0] = fspc_0
            fspc[set0] = fspc_0

        # ----------------------MPC recovery (n -> g)----------------------
        if GMN is not None:
            # Expand n-set solution to g-set: xg = GMN @ xn
            xn_solve = np.where(np.isnan(xn_full), 0.0, xn_full)
            xg = np.asarray(GMN @ xn_solve).ravel()
            # SPC forces in g-set: fspc_g = Kgg_orig @ xg - Fg_orig
            fspc_g = np.full(ndof, 0.0, dtype=fdtype)
            # Rebuild Fg in g-set
            Fg_g = np.asarray(Kgg_orig @ xg).ravel()
            # For the SPC DOFs in the original g-set
            sset_orig, sset_b_orig, unused_xg_orig2 = _build_xg(
                model, dof_map, ndof, subcase)
            fspc_g[sset_orig] = Fg_g[sset_orig]
            xg_out = xg
            Fg_out = Fg_g
            fspc = fspc_g
            ndof_out = ndof

            # MPC forces: F_mpc = Kgg @ xg - Fg_applied at dependent DOFs
            # The MPC force is the reaction needed to enforce the constraint
            Fg_applied = np.where(np.isnan(Fg_gset), 0.0, Fg_gset)
            F_internal = np.asarray(Kgg_orig @ xg).ravel()
            fmpc = np.zeros(ndof, dtype=fdtype)
            fmpc[mset] = F_internal[mset] - Fg_applied[mset]
            self.fmpc = fmpc
            log.debug(f"  MPC forces at {len(mset)} dependent DOFs")
        else:
            xg_out = xn_full
            Fg_out = Fg
            ndof_out = ndof
            self.fmpc = np.zeros(ndof, dtype=fdtype)

        log.debug(f"xa = {xa}")
        log.debug(f"xs = {xs}")
        log.debug(f"fspc = {fspc}")

        xg = xg_out
        Fg = Fg_out
        self.xg = xg
        self.Fg = Fg

        xb = xg_to_xb(model, xg, ngrid, ndof_per_grid)
        Fb = xg_to_xb(model, Fg, ngrid, ndof_per_grid)

        # log.debug(f'Fs = {Fs}')
        # log.debug(f'Fb = {Fb}')
        # log.debug(f'xb = {xb}')

        self._save_displacment(
            f06_file,
            subcase,
            itime,
            ntimes,
            node_gridtype,
            xg,
            ngrid,
            ndof_per_grid,
            title=title,
            subtitle=subtitle,
            label=label,
            fdtype=fdtype,
            page_num=page_num,
            page_stamp=page_stamp,
        )

        self._save_applied_load(
            f06_file,
            subcase,
            itime,
            ntimes,
            node_gridtype,
            Fg_oload,
            ngrid,
            ndof_per_grid,
            title=title,
            subtitle=subtitle,
            label=label,
            fdtype=fdtype,
            page_num=page_num,
            page_stamp=page_stamp,
        )

        self._save_spc_forces(
            f06_file,
            subcase,
            itime,
            ntimes,
            node_gridtype,
            fspc,
            ngrid,
            ndof_per_grid,
            title=title,
            subtitle=subtitle,
            label=label,
            fdtype=fdtype,
            page_num=page_num,
            page_stamp=page_stamp,
        )

        self._save_mpc_forces(
            f06_file,
            subcase,
            itime,
            ntimes,
            node_gridtype,
            self.fmpc,
            ngrid,
            ndof_per_grid,
            title=title,
            subtitle=subtitle,
            label=label,
            fdtype=fdtype,
            page_num=page_num,
            page_stamp=page_stamp,
        )

        # SPCFORCE resultant
        if 'SPCFORCES' in subcase:
            _write_spcforce_resultant(
                f06_file, fspc, ngrid, ndof_per_grid, isubcase,
                page_stamp, page_num, log)

        if 'GPFORCE' in subcase:
            _write_gpforce_balance(
                f06_file, self.op2, node_gridtype, Fg_oload, fspc,
                Kgg_orig, xg, isubcase, ngrid, ndof_per_grid,
                title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp, page_num=page_num)

        op2 = self.op2
        page_stamp += "\n"
        recover_force_101(
            f06_file,
            op2,
            self.model,
            dof_map,
            subcase,
            xb,
            title=title,
            subtitle=subtitle,
            label=label,
            page_stamp=page_stamp,
        )

        recover_strain_101(
            f06_file, op2, self.model, dof_map, subcase, xb,
            title=title, subtitle=subtitle, label=label,
            page_stamp=page_stamp,
        )
        recover_stress_101(
            f06_file, op2, self.model, dof_map, subcase, xb,
            title=title, subtitle=subtitle, label=label,
            page_stamp=page_stamp,
        )
        recover_strain_energy_101(
            f06_file,
            op2,
            self.model,
            dof_map,
            subcase,
            xb,
            title=title,
            subtitle=subtitle,
            label=label,
            page_stamp=page_stamp,
        )
        self.log.info("finished")
        out = {}
        return out, page_num, end_options

    def _save_displacment(
        self,
        f06_file: TextIO,
        subcase: Subcase,
        itime: int,
        ntimes: int,
        node_gridtype: NDArrayN2int,
        xg: NDArrayNfloat,
        ngrid: int,
        ndof_per_grid: int,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        fdtype: str = "float32",
        page_num: int = 1,
        page_stamp: str = "PAGE %s",
    ) -> int:
        f06_request_name = "DISPLACEMENT"
        table_name = "OUGV1"
        # self.log.debug(f'xg = {xg}')
        page_num = save_static_table(
            f06_file,
            subcase,
            itime,
            ntimes,
            node_gridtype,
            xg,
            RealDisplacementArray,
            f06_request_name,
            table_name,
            self.op2.displacements,
            ngrid,
            ndof_per_grid,
            title=title,
            subtitle=subtitle,
            label=label,
            fdtype=fdtype,
            page_num=page_num,
            page_stamp=page_stamp,
        )
        return page_num

    def _save_spc_forces(
        self,
        f06_file: TextIO,
        subcase: Subcase,
        itime: int,
        ntimes: int,
        node_gridtype: NDArrayN2int,
        fspc: NDArrayNfloat,
        ngrid: int,
        ndof_per_grid: int,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        fdtype: str = "float32",
        page_num: int = 1,
        page_stamp: str = "PAGE %s",
    ) -> int:
        f06_request_name = "SPCFORCES"
        table_name = "OQG1"
        # self.log.debug(f'Fg = {Fg}')
        page_num = save_static_table(
            f06_file,
            subcase,
            itime,
            ntimes,
            node_gridtype,
            fspc,
            RealSPCForcesArray,
            f06_request_name,
            table_name,
            self.op2.spc_forces,
            ngrid,
            ndof_per_grid,
            title=title,
            subtitle=subtitle,
            label=label,
            fdtype=fdtype,
            page_num=page_num,
            page_stamp=page_stamp,
        )
        return page_num

    def _save_mpc_forces(
        self,
        f06_file: TextIO,
        subcase: Subcase,
        itime: int,
        ntimes: int,
        node_gridtype: NDArrayN2int,
        fmpc: NDArrayNfloat,
        ngrid: int,
        ndof_per_grid: int,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        fdtype: str = "float32",
        page_num: int = 1,
        page_stamp: str = "PAGE %s",
    ) -> int:
        """Save MPC forces to F06 and OP2."""
        from pyNastran.op2.tables.oqg_constraintForces.oqg_mpc_forces import (
            RealMPCForcesArray,
        )
        from .recover.utils import get_plot_request
        from .utils import recast_data

        f06_request_name = "MPCFORCES"
        unused_nids_write, write_f06, write_op2, quick_return = get_plot_request(
            subcase, f06_request_name)
        if quick_return:
            return page_num

        idtype2, fdtype2 = recast_data("int32", fdtype)
        isubcase = subcase.id
        nnodes = node_gridtype.shape[0]
        data = np.zeros((ntimes, nnodes, 6), dtype=fdtype2)
        ngrid_dofs = ngrid * ndof_per_grid
        _fgi = fmpc[:ngrid_dofs].reshape(ngrid, ndof_per_grid)
        data[itime, :ngrid, :] = _fgi
        data[itime, ngrid:, 0] = fmpc[ngrid_dofs:]

        # Use OQG1 for data_code creation, then override table_name for write_f06
        table_name = "OQG1"
        mpc_obj = RealMPCForcesArray.add_static_case(
            table_name, node_gridtype, data, isubcase,
            is_sort1=True, is_random=False, is_msc=True,
            random_code=0, title=title, subtitle=subtitle, label=label)
        mpc_obj.table_name = "OQMG1"

        if write_f06:
            page_num = mpc_obj.write_f06(
                f06_file, header=None,
                page_stamp=page_stamp, page_num=page_num,
                is_mag_phase=False, is_sort1=True)
            f06_file.write("\n")
        if write_op2:
            self.op2.mpc_forces[isubcase] = mpc_obj
        return page_num

    def _save_applied_load(
        self,
        f06_file: TextIO,
        subcase: Subcase,
        itime: int,
        ntimes: int,
        node_gridtype: NDArrayN2int,
        Fg: NDArrayNfloat,
        ngrid: int,
        ndof_per_grid: int,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        fdtype: str = "float32",
        page_num: int = 1,
        page_stamp: str = "PAGE %s",
    ) -> int:
        f06_request_name = "OLOAD"
        table_name = "OPG1"
        # self.log.debug(f'Fg = {Fg}')
        page_num = save_static_table(
            f06_file,
            subcase,
            itime,
            ntimes,
            node_gridtype,
            Fg,
            RealLoadVectorArray,
            f06_request_name,
            table_name,
            self.op2.load_vectors,
            ngrid,
            ndof_per_grid,
            title=title,
            subtitle=subtitle,
            label=label,
            fdtype=fdtype,
            page_num=page_num,
            page_stamp=page_stamp,
        )
        return page_num

    def run_sol_103_modes(
        self,
        subcase: Subcase,
        f06_file: TextIO,
        page_stamp: str,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        page_num: int = 1,
        idtype: str = "int32",
        fdtype: str = "float64",
        write_op2: bool = True,
    ):
        """
        [M]{xdd} + [C]{xd} + [K]{x} = {F}
        [M]{xdd} + [K]{x} = {0}
        -[M]{xdd}λ^2 + [K]{x} = {0}
        {X}(λ^2 - [M]^-1[K]) = {0}
        λ^2 - [M]^-1[K] = {0}
        λ^2 = [M]^-1[K]
        [A][X] = [X]λ^2
        """
        model = self.model
        log = model.log
        log.debug(f"run_sol_103 (modes)")
        assert len(model.methods), "SOL 103 (modes) requires a METHOD and a EIGR/EIGRL card"
        end_options = [
            "SEMR",  # MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES)
            "SEKR",  # STIFFNESS MATRIX REDUCTION STEP
            "MODES",  # run modes
        ]
        op2 = self.op2
        # write_f06 = True

        node_gridtype = _get_node_gridtype(model, idtype=idtype)

        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof = get_ndof(self.model, subcase)
        # Build GMN for MPC reduction
        GMN, mset = self.build_GMN(subcase, dof_map, ndof, fdtype=fdtype)
        self.GMN = GMN
        self.mset = mset

        out = _run_modes(
            model,
            subcase,
            ngrid,
            ndof_per_grid,
            ndof,
            node_gridtype,
            dof_map,
            Kgg_override=self.Kgg_override,
            GMN=GMN,
            mset=mset,
            idtype=idtype,
            fdtype=fdtype,
        )
        phig = out["modes_phig"]
        eigenvalue = out["modes_eigenvalue"]
        Mhh = out["modes_Mhh"]
        Khh = out["modes_Khh"]
        # out['modes_phigg'] = xg_out
        nmode = len(eigenvalue)

        isubcase = subcase.id
        mode_cycle = eigenvalue
        eigenvalue_obj = RealEigenvalues(title, "LAMA", nmodes=nmode)

        cycle = np.sqrt(np.abs(eigenvalue)) / (2.0 * np.pi)
        radian = np.sqrt(np.abs(eigenvalue))

        eigenvalue_obj.mode = np.arange(nmode, dtype="int32") + 1
        eigenvalue_obj.extraction_order = np.arange(nmode, dtype="int32") + 1
        eigenvalue_obj.eigenvalues = eigenvalue
        eigenvalue_obj.radians = radian
        eigenvalue_obj.cycles = cycle
        eigenvalue_obj.generalized_mass = np.diag(Mhh)
        eigenvalue_obj.generalized_stiffness = np.diag(Khh)
        op2.eigenvalues[title] = eigenvalue_obj
        # op2.eigenvalues[isubcase] = eigenvalues_obj

        # Mass participation factors
        Mgg = out["Mgg"]
        nnode_g = node_gridtype.shape[0]
        mpf = compute_mass_participation(phig, Mgg, nnode_g, nmode)
        out["mass_participation"] = mpf

        _write_mass_participation_f06(f06_file, mpf, cycle, nmode)
        log.info(f"mass participation cumulative (Tx): {mpf['cumulative_ratio'][-1, 0]:.4f}")

        (
            write_phi_f06,
            write_phi_op2,
            write_phi_pch,
            unused_options,
            phi_set,
        ) = get_f06_op2_pch_set(subcase, "DISPLACEMENT")
        write_eigenvector = any((write_phi_f06, write_phi_op2))

        modes = np.arange(1, nmode + 1, dtype=idtype)
        eigenvalue = eigenvalue.astype("float32")

        if write_eigenvector:
            nnode = node_gridtype.shape[0]
            # (1,2050,18)
            # print(xg_out.shape)
            node_gridtypei, phi, nnodei = slice_modal_set(
                node_gridtype, phig, nnode, nmode, phi_set
            )

            data = phi.reshape((nmode, nnodei, 6)).astype("float32")
            table_name = "OUGV1"
            eigenvector_obj = RealEigenvectorArray.add_modal_case(
                table_name,
                node_gridtypei,
                data,
                isubcase,
                modes,
                eigenvalue,
                mode_cycle,
                is_sort1=True,
                is_random=False,
                is_msc=True,
                random_code=0,
                title=title,
                subtitle=subtitle,
                label=label,
            )
            eigenvector_obj.nonlinear_factor = 1
            op2.eigenvectors[isubcase] = eigenvector_obj

        write_eigenvalue_f06 = True
        write_eigenvalue_op2 = True
        if write_eigenvalue_f06:
            str(page_num)
            header = []
            page_num = eigenvalue_obj.write_f06(
                f06_file, header=header, page_stamp=page_stamp, page_num=page_num
            )
            f06_file.write("\n")

        if write_phi_f06:
            page_num = eigenvector_obj.write_f06(
                f06_file,
                header=None,
                page_stamp=page_stamp,
                page_num=page_num,
                is_mag_phase=False,
                is_sort1=True,
            )
            f06_file.write("\n")
        # fspc = Ksa @ xa + Kss @ xs
        # Fs[ipositive] = Fsi

        # Fg[aset] = Fa
        # Fg[sset] = fspc
        # log.debug(f'xa = {xa}')
        # log.debug(f'Fs = {Fs}')
        # log.debug(f'xg = {xg}')
        # log.debug(f'Fg = {Fg}')

        op2 = self.op2
        page_stamp += "\n"

        # TODO: add transform (need for rods, quads, etc., but not springs)
        phib = phig
        recover_force_103(
            f06_file,
            op2,
            self.model,
            dof_map,
            subcase,
            phig,
            phib,
            modes,
            eigenvalue,  # freqs,
            title=title,
            subtitle=subtitle,
            label=label,
            page_stamp=page_stamp,
        )

        # Per-mode force/stress/strain/ESE recovery using static routines
        nmode = phig.shape[0]
        for imode in range(nmode):
            xb_mode = phib[imode, :]
            recover_force_101(
                f06_file, op2, self.model, dof_map, subcase,
                xb_mode, title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp,
            )
            recover_strain_101(
                f06_file, op2, self.model, dof_map, subcase,
                xb_mode, title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp,
            )
            recover_stress_101(
                f06_file, op2, self.model, dof_map, subcase,
                xb_mode, title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp,
            )
            recover_strain_energy_101(
                f06_file, op2, self.model, dof_map, subcase,
                xb_mode, title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp,
            )
        str(f06_file)
        str(page_stamp)
        return out, page_num, end_options

    def run_sol_111_modal_frequency(
        self,
        subcase: Subcase,
        f06_file: TextIO,
        page_stamp: str,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        page_num: int = 1,
        idtype: str = "int32",
        fdtype: str = "float64",
    ):
        model = self.model
        log = model.log
        isubcase = subcase.id
        op2 = self.op2
        # -----------------------------------------------------------------------
        node_gridtype = _get_node_gridtype(model, idtype=idtype)

        # handles MAX/MASS normalization
        out, page_num, end_options = self.run_sol_103_modes(
            subcase,
            f06_file,
            page_stamp=page_stamp,
            title=title,
            subtitle=subtitle,
            label=label,
            page_num=page_num,
            idtype=idtype,
            fdtype=fdtype,
        )
        end_options.append("DYNRED")
        # aset = out['aset']
        Kgg = out["Kgg"]
        # Mgg = out['Mgg']
        eigenvalue = out["modes_eigenvalue"]
        # print(list(out))
        phit = out["modes_phig"]  # phi-gg

        nmode = len(eigenvalue)
        ndof_g = Kgg.shape[0]
        # ndof_a = Kaa.shape[0]

        # natural frequency
        omegan = np.sqrt(np.abs(eigenvalue))

        # frequency
        freq = get_frequencies(model, subcase, omegan)
        nfreq = len(freq)
        omega = 2 * np.pi * freq

        # TODO: parse inputs for loads
        Fg = np.ones((ndof_g, 1), dtype="complex64")

        Cstr_gg, zomegan2 = get_freq_damping(model, omegan, ndof_g)

        # nfreq = len(omegas)
        # shape = (ndof, ndof, nfreq)

        # phiT = phi.T
        #  modal space (h); sometimes called x
        # print(phi.shape, Mgg.shape)
        nmodes, norm_str = get_real_eigenvalue_method(model, subcase)
        # phi = phit.T
        Mhh = out["modes_Mhh"]
        Khh = out["modes_Khh"]
        if norm_str == "MAX":
            raise RuntimeError("norm_str=MAX and should be MASS (it makes the math harder)")
        # Mhh = phit @ Mgg @ phi
        # Chh = phi @ Cgg @ phi.t
        # Khh = phi @ Kgg @ phi.T
        # Fh = phi @ Fg

        Fh = phit @ Fg
        # print('Fh:\n', Fh)
        omega2 = omega * omega  # frequncy
        omegan2 = np.diag(Khh)  # natural frequency squared
        omegan = np.sqrt(omegan2)

        k = np.diag(Khh)  # TODO: this part is weird
        p = 1.0
        # print(f'ndof_g={ndof_g} phi.shape={phi.shape}')
        xq = np.zeros((nfreq, nmode), dtype="complex64")
        for ifreq, omega2i in enumerate(omega2):
            log.info(f"ifreq={ifreq}")
            # denom = np.sqrt((1-omega2i/omegan2)**2 + zomegan2**2)
            # p / (k * denom) * np.sin(omegan)

            # $$ u(t) =
            tf = np.sqrt((1 - omega2i / omegan2) ** 2 + (zomegan2 / omegan) ** 2)
            mag = p / k * tf
            phase = -np.atan2(zomegan2 / omegan, 1 - omega2i / omegan2)
            real = mag * np.cos(phase)
            imag = mag * np.sin(phase)
            xq_freq = real + 1j * imag
            # print(xq_freq.shape)
            xq[ifreq, :] = xq_freq
            log.info(xq_freq)

        xg = xq @ phit  # (1000,2) * (2,18) = (1000,18)
        # vg = (1j * omega)[:, np.newaxis] * xg[np.newaxis, :]
        # ValueError: operands could not be broadcast together with shapes (2050,) (2050,18)
        # print(omega2.shape, xg.shape)
        # ag = -omega2[:, np.newaxis] * xg[np.newaxis, :]

        key_to_factor = {
            "DISPLACEMENT": (op2.displacements, ComplexDisplacementArray, False, 1.0),
            "VELOCITY": (op2.velocities, ComplexVelocityArray, True, 1j * omega),
            "ACCELERATION": (op2.accelerations, ComplexAccelerationArray, True, -omega * omega),
        }
        header = ["\n", "\n", "\n"]
        nnode = node_gridtype.shape[0]
        for key, (slot, obj_class, apply_scale, factor) in key_to_factor.items():
            write_f06, write_op2, write_pch, options, set_data = get_f06_op2_pch_set(subcase, key)
            write_data = np.any((write_f06, write_op2))
            if not write_data:
                log.warning(f"skipping {key!r}")
                continue
            if not apply_scale:  # displacment
                data = xg.copy()
                assert data.shape == (nfreq, ndof_g), data.shape
            else:  # velocity, acceleration
                data = factor[:, np.newaxis] * xg[np.newaxis, :]
                data = data.reshape(nfreq, ndof_g)
                assert data.shape == (nfreq, ndof_g), (data.shape, (nfreq, ndof_g))

            # print(key, data.shape, (nfreq, nnode))
            # (2050,18)
            node_gridtypei, data, nnodei = slice_freq_set(
                node_gridtype, data, nnode, nfreq, set_data
            )
            assert data.shape == (nfreq, nnodei, 6), data.shape

            # eigenvaluef = eigenvalue.astype('float32')
            table_name = "OUGV1"
            obj = obj_class.add_freq_case(
                table_name,
                node_gridtypei,
                data,
                isubcase,
                freq,
                is_sort1=True,
                is_random=False,
                is_msc=True,
                random_code=0,
                title=title,
                subtitle=subtitle,
                label=label,
            )
            # obj.nonlinear_factor = 1
            slot[isubcase] = obj

            header = ["", "", ""]
            if write_f06:
                is_mag_phase = get_mag_phase_from_options(options)
                page_num = obj.write_f06(
                    f06_file,
                    header=header,
                    page_stamp=page_stamp,
                    page_num=page_num,
                    is_mag_phase=is_mag_phase,
                    is_sort1=True,
                )
                f06_file.write("\n")

        # -w^2 [Mhat]qdd + j*w[Bhat]qd + [Khat]q = F

        # ------------------
        # SPC forces ???
        # -w^2 [Mhat]qdd + j*w[Bhat]qd + [Khat]q = F
        # -w^2 [Mhat]q + [Khat]q = 0
        # TODO: how do I get the modal constraint force?

        # ------------------
        dof_map = None
        recover_force_freq(
            f06_file,
            op2,
            self.model,
            dof_map,
            subcase,
            xq,
            # modes, eigenvalue,
            freq,
            title=title,
            subtitle=subtitle,
            label=label,
            page_stamp=page_stamp,
        )
        return out, page_num, end_options

    def run_sol_105_buckling(
        self,
        subcase: Subcase,
        f06_file: TextIO,
        page_stamp: str,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        page_num: int = 1,
        idtype: str = "int32",
        fdtype: str = "float64",
        write_op2: bool = True,
    ):
        """SOL 105 linear buckling: (K + lambda*KD)*x = 0.

        Subcase must have LOAD, SPC, and METHOD.
        STATSUB(PRELOAD) is not yet supported — uses the same subcase for preload.
        """
        model = self.model
        log = model.log
        log.debug("run_sol_105 (buckling)")

        end_options = ["SEKR", "MODES"]

        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof = get_ndof(model, subcase)

        from .build_stiffness import build_Kgg

        if self.Kgg_override is not None:
            Kgg_arr = self.Kgg_override
            if issparse(Kgg_arr):
                Kgg = Kgg_arr
            else:
                Kgg = csc_matrix(Kgg_arr[:ndof, :ndof])
        else:
            Kgg = build_Kgg(model, dof_map, ndof, ngrid, ndof_per_grid)

        # Static preload solve
        sset, sset_b, xg = _build_xg(model, dof_map, ndof, subcase)
        Fb = self.build_Fb(xg, sset_b, dof_map, ndof, subcase)

        free_dofs = np.where(~sset_b)[0]
        Kff = Kgg.tocsc()[np.ix_(free_dofs, free_dofs)].toarray()
        Ff = Fb[free_dofs]

        u_free = np.linalg.solve(Kff, Ff)
        u_global = np.zeros(ndof, dtype=fdtype)
        u_global[free_dofs] = u_free

        log.debug("  static preload solved")

        # Build geometric stiffness from preload stress state
        from .shells import build_KDgg_cquad4, build_KDgg_ctria3
        from .build_stiffness import build_KDgg_beam
        from .solids import build_KDgg_solids

        KDgg = dok_matrix((ndof, ndof))
        build_KDgg_cquad4(model, KDgg, dof_map, u_global)
        build_KDgg_ctria3(model, KDgg, dof_map, u_global)
        build_KDgg_beam(model, KDgg, dof_map, u_global)
        build_KDgg_solids(model, KDgg, dof_map, u_global)

        KDff = KDgg.tocsc()[np.ix_(free_dofs, free_dofs)].toarray()

        log.debug("  geometric stiffness assembled")

        # Solve buckling eigenproblem: (K + lambda*KD)*x = 0
        # => K*x = -lambda*KD*x
        # Use scipy generalized eigenvalue: Kff @ x = lambda * (-KDff) @ x
        from scipy.linalg import eigh

        neigenvalue, _ = get_real_eigenvalue_method(model, subcase)
        neg_KDff = -KDff

        # Use eigh for symmetric positive-definite B (standard buckling)
        # If -KD isn't positive definite, fall back to general eig
        try:
            eigenvalues_all, eigvecs_all = eigh(Kff, neg_KDff)
        except np.linalg.LinAlgError:
            from scipy.linalg import eig
            eigenvalues_all, eigvecs_all = eig(Kff, neg_KDff)
            eigenvalues_all = eigenvalues_all.real
            eigvecs_all = eigvecs_all.real

        # Keep only positive eigenvalues (physical buckling modes) sorted ascending
        pos_mask = eigenvalues_all > 0
        eigenvalues_pos = eigenvalues_all[pos_mask]
        eigvecs_pos = eigvecs_all[:, pos_mask]
        sort_idx = np.argsort(eigenvalues_pos)
        eigenvalue = eigenvalues_pos[sort_idx[:neigenvalue]]
        xa_ = eigvecs_pos[:, sort_idx[:neigenvalue]]
        nmode = len(eigenvalue)

        # Buckling eigenvalues (critical load factors)
        nmode = len(eigenvalue)
        log.info(f"  buckling eigenvalues (critical load factors): {eigenvalue[:min(5,nmode)]}")

        # Store results
        log.info(f"  buckling eigenvalues: {eigenvalue[:min(5, nmode)]}")

        out = {
            "buckling_eigenvalues": eigenvalue,
            "buckling_modes": xa_,
            "free_dofs": free_dofs,
        }
        self.buckling_eigenvalues = eigenvalue
        return out, page_num, end_options

    def run_sol_31_craig_bampton(
        self,
        subcase: Subcase,
        f06_file: TextIO,
        page_stamp: str,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        page_num: int = 1,
        idtype: str = "int32",
        fdtype: str = "float64",
    ):
        """SOL 31 Craig-Bampton (GEN CB MODEL).

        Boundary DOFs defined by SUPORT/SUPORT1 cards.
        Number of modes from METHOD (EIGRL/EIGR).
        """
        from .craig_bampton import run_craig_bampton

        model = self.model
        log = model.log
        log.debug("run_sol_31 (Craig-Bampton)")

        end_options = ["SEMG", "SEMR", "SEKR"]

        model.setup(run_geom_check=True)
        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof = get_ndof(model, subcase)

        # Build global matrices
        if self.Kgg_override is not None:
            from scipy.sparse import issparse as _issparse
            Kgg_in = self.Kgg_override
            if _issparse(Kgg_in):
                Kgg = Kgg_in[:ndof, :ndof].tocsc()
            else:
                Kgg = csc_matrix(Kgg_in[:ndof, :ndof])
        else:
            Kgg = build_Kgg(model, dof_map, ndof, ngrid, ndof_per_grid,
                            idtype=idtype, fdtype=fdtype)

        Mbb = build_Mbb(model, subcase, dof_map, ndof, fdtype=fdtype)
        from .build_stiffness import Kbb_to_Kgg
        Mgg = Kbb_to_Kgg(model, Mbb, ngrid, ndof_per_grid, inplace=False)

        # Apply SPC constraints — reduce to free (f) set
        gset_b = ps_to_sg_set(ndof, ps)
        sset, sset_b, xg = _build_xg(model, dof_map, ndof, subcase)
        free_dofs = np.where(~sset_b)[0]

        Kff = Kgg.tocsc()[np.ix_(free_dofs, free_dofs)]
        Mff = Mgg.tocsc()[np.ix_(free_dofs, free_dofs)]

        # Build DOF map for the free set
        # Re-map dof_map indices to free-set indices
        g_to_f = np.full(ndof, -1, dtype="int32")
        g_to_f[free_dofs] = np.arange(len(free_dofs), dtype="int32")

        # Identify R-set DOFs from SUPORT/SUPORT1 cards
        r_set_dofs = []
        suport = model.suport
        if suport.n > 0:
            for i in range(suport.n):
                nid = int(suport.node_id[i])
                comp_str = str(suport.component[i])
                for c in comp_str:
                    dof = int(c)
                    g_idx = dof_map[(nid, dof)]
                    f_idx = g_to_f[g_idx]
                    if f_idx >= 0:
                        r_set_dofs.append((nid, dof))

        suport1 = model.suport1
        if suport1.n > 0:
            for i in range(suport1.n):
                nid = int(suport1.node_id[i])
                comp_str = str(suport1.component[i])
                for c in comp_str:
                    dof = int(c)
                    g_idx = dof_map[(nid, dof)]
                    f_idx = g_to_f[g_idx]
                    if f_idx >= 0:
                        r_set_dofs.append((nid, dof))

        if not r_set_dofs:
            raise RuntimeError(
                "SOL 31 requires SUPORT or SUPORT1 cards to define boundary DOFs"
            )

        # Build DOF map for the free set
        f_dof_map: DOF_MAP = {}
        for (nid, dof), g_idx in dof_map.items():
            f_idx = g_to_f[g_idx]
            if f_idx >= 0:
                f_dof_map[(nid, dof)] = f_idx

        # Number of eigenvalues
        neigenvalues = 20
        if "METHOD" in subcase:
            neigenvalues, _ = get_real_eigenvalue_method(model, subcase)

        log.info(f"  R-set DOFs: {len(r_set_dofs)}, modes requested: {neigenvalues}")

        # Run Craig-Bampton
        cb_result = run_craig_bampton(
            Kff, Mff, f_dof_map, r_set_dofs,
            neigenvalues=neigenvalues, log=log,
        )

        # Write summary to F06
        eigenvalues = cb_result["eigenvalues"]
        nvec = len(eigenvalues)
        f06_file.write(f"\n{'='*72}\n")
        f06_file.write(f"  CRAIG-BAMPTON MODEL GENERATION (SOL 31)\n")
        f06_file.write(f"{'='*72}\n\n")
        f06_file.write(f"  Number of boundary (R-set) DOFs: {len(r_set_dofs)}\n")
        f06_file.write(f"  Number of fixed-interface modes: {nvec}\n")
        f06_file.write(f"  Total CB DOFs: {cb_result['num_cb_dofs']}\n\n")

        f06_file.write(f"  FIXED-INTERFACE NATURAL FREQUENCIES\n")
        f06_file.write(f"  {'MODE':>4s}  {'EIGENVALUE':>14s}  {'FREQUENCY (Hz)':>14s}\n")
        for i, ev in enumerate(eigenvalues):
            freq_hz = np.sqrt(abs(ev)) / (2.0 * np.pi)
            f06_file.write(f"  {i+1:4d}  {ev:14.6e}  {freq_hz:14.6f}\n")
        f06_file.write("\n")

        # Store results on solver object
        self.cb_result = cb_result

        # Write CB matrices to OP4 and HDF5
        from .craig_bampton import write_cb_to_op4, write_cb_to_h5
        base_name = os.path.splitext(model.bdf_filename)[0]
        op4_filename = base_name + ".cb.op4"
        h5_filename = base_name + ".cb.h5"
        write_cb_to_op4(op4_filename, cb_result, is_binary=True, precision="double")
        log.info(f"  CB matrices written to: {op4_filename}")
        try:
            write_cb_to_h5(h5_filename, cb_result)
            log.info(f"  CB matrices written to: {h5_filename}")
        except ImportError:
            log.warning("  pytables not installed; skipping HDF5 output")

        out = cb_result
        return out, page_num, end_options

    def run_sol_108_direct_frequency(
        self,
        subcase: Subcase,
        f06_file: TextIO,
        page_stamp: str,
        title: str = "",
        subtitle: str = "",
        label: str = "",
        page_num: int = 1,
        idtype: str = "int32",
        fdtype: str = "float64",
    ):
        """
        Direct frequency response

        frequency
        [M]{xdd} + [C]{xd} + [K]{x} = {F}
        {φ} ([M]s^2 + [C]{s} + [K]) = {F}
        {φ} ([M]s^2 + [C]{s} + [K]) = {F}
        {φ} ([M]s^2 + [K]) = {F}

        {d}*sin(wt) ([M]s^2 + [K]) = {0}
        [M]{xdd} + [K]{x} = {0}
        x = A*sin(wt)
        xd = dx/dt = A*w*cos(wt)
        xdd = -A*w^2 * sin(wt)
        -w^2[M]{A}*sin(wt) + [K]{A}*sin(wt) = {0}
        (-w^2[M] + [K]){A} = {0}
        -w^2[M] + [K] = 0
        w^2 = [M]^-1 * [K]

        -w^2[M]{A}*sin(wt) + [K]{A}*sin(wt) = {F}sin(wt)
        (-w^2[M] + [K]){A} = {F}
        {A} = (-w^2[M] + [K])^-1 * {F}
        {A} = ([K] - w^2[M])^-1 * {F}

        {φ}^T[M]{φ}{qdd} + {φ}^T[B]{φ}{qd} {φ}^T[K]{φ}{q} = {φ}^T{P(t)}
        {qdd} + 2*zeta*w_n{qd} w_n^2{q} = 1/M_n {φ}^T{P(t)}
        {u_total} = {u_dynamic} + {u_static} = [PHI]{q} + {u_static}

        """
        str(subcase)
        str(f06_file)
        str(fdtype)
        str(idtype)
        # https://pyfrf.readthedocs.io/en/latest/Showcase.html
        # ndof = K.shape[0]

        # end_options = [
        #'SEMR',  # MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES)
        #'SEKR',  # STIFFNESS MATRIX REDUCTION STEP
        # ]
        model = self.model
        # op2 = self.op2
        # ---------------------------------------------------
        # title = ''
        # today = None
        # page_stamp = self.op2.make_stamp(title, today) # + '\n'

        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof = get_ndof(self.model, subcase)

        GMN, mset = self.build_GMN(subcase, dof_map, ndof, fdtype=fdtype)
        self.GMN = GMN
        self.mset = mset

        out = _run_modes(
            model,
            subcase,
            ngrid,
            ndof_per_grid,
            ndof,
            node_gridtype,
            dof_map,
            GMN=GMN,
            mset=mset,
            idtype=idtype,
            fdtype=fdtype,
        )

        aset = out["aset"]
        Kaa = out["Kaa"]
        Maa = out["Maa"]
        eigenvalue = out["modes_eigenvalue"]
        omegas = np.sqrt(np.abs(eigenvalue))

        freqs = get_frequencies(model, subcase, omegas)
        nfreq = len(freqs)
        omegai = 2 * np.pi * freqs

        ndof_a = Kaa.shape[0]
        Caa = np.eye(ndof_a)
        np.fill_diagonal(Caa, 0.01)

        # nfreq = len(omegas)
        shape = (ndof, ndof, nfreq)
        tf_matrix = np.zeros(shape, dtype="complex128")  # full system 3x3 TF matrix
        tf_mask = np.ones((ndof, ndof, nfreq), dtype="bool")
        tf_mask[~aset, :, :] = 0
        tf_mask[:, ~aset, :] = 0
        tf_matrix_a = np.ones((ndof_a, ndof_a, nfreq), dtype="complex128")
        # tf_matrix_a = np.ma.masked_array(tf_matrix, tf_mask)
        backend = get_solver()
        for i, omegai in enumerate(omegas):
            mat = -(omegai**2) * Maa + 1j * omegai * Caa + Kaa
            res = backend.inv(mat)
            tf_matrix_a[:, :, i] = res

        tf_matrix = np.zeros(shape, dtype="complex128")  # full system 3x3 TF matrix
        tf_matrix.ravel()[tf_mask.ravel()] = tf_matrix_a.ravel()

        # TODO: save tf_matrix to ACCEL response
        # return out, tf_matrix
        end_options = [
            "SEMG",  # STIFFNESS AND MASS MATRIX GENERATION STEP
            "SEMR",  # MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES)
            "SEKR",  # STIFFNESS MATRIX REDUCTION STEP
            "SELG",  # LOAD MATRIX GENERATION STEP
            "SELR",  # LOAD MATRIX REDUCTION STEP
        ]
        return out, page_num, end_options


def remove_rows(Kgg: NDArrayNNfloat, aset: NDArrayNint, idtype: str = "int32") -> NDArrayNNfloat:
    """
    Applies AUTOSPC to the model (sz)

    mp  DOFs eliminated by multipoint constraints.
    mr  DOFs eliminated by multipoint constraints created by the rigid
        elements using the LGELIM method on the Case Control command RIGID.
    sb* DOFs eliminated by single-point constraints that are included
        in boundary condition changes and by the AUTOSPC feature.
        (See the sz set)
    sg* DOFs eliminated by single-point constraints that are specified
        on the PS field on GRID Bulk Data entries.
    sz  DOFs eliminated by the AUTOSPC feature.
    o   DOFs omitted by structural matrix partitioning.
    q   Generalized DOFs assigned to component modes and residual vectors.
    r   Reference DOFs used to determine free body motion.
    c   DOFs that are free during component mode synthesis or dynamic reduction.
    b   DOFs fixed during component mode analysis or dynamic reduction.
    lm  Lagrange multiplier DOFs created by the rigid elements
        using the LAGR method on the Case Control command, RIGID.
    e   Extra DOFs introduced in dynamic analysis.
    sa  Permanently constrained aerodynamic DOFs.
    k   Aerodynamic mesh point set for forces and displacements on the aero mesh.
    j   Aerodynamic mesh collocation point set (exact physical
        interpretation is dependent on the aerodynamic theory).


    Input Sets
    ----------
    mp  DOFs eliminated by multipoint constraints.
    mr  DOFs eliminated by multipoint constraints created by the rigid
        elements using the LGELIM method on the Case Control command RIGID.
    g   all DOFs including scalar DOFs
    o   Omit set
    s   SPC sets (sb + sg)

    Source of Sets
    --------------
    g     6-DOFs per GRID, 1-DOF SPOINTs,
          Nmodes-DOF EPOINTs (Extra Dynamic/Modal Point; ),
          nRBE-DOF LPOINTs ("Lagrange" Point; e.g., RBE3),
          nHarmonic-DOF HPOINTs (Harmonic Point; e.g., RINGFL)
    m     MPC, MPCADD, RBAR/RBAR1, RBE1/RBE2/RBE3, RROD, RSPLINE,
          RJOINT, RSSCON, RTRPLT/RTRPLT1, GMBC, GMSPC*
    r     SUPORT/SUPORT1, SUPAX
    s     SPC/SPC1, SPCADD, BNDGRID, AUTOSPC
    sb    SPC/SPC1, SPCADD, SPCAX, FLSYM, GMSPC*, BNDGRID, PARAM,AUTOSPC,YES
    sz    PARAM,AUTOSPC,YES
    sg    GRID, GRIDB, GRDSET (PS field)
    o     OMIT/OMIT1, GRID (SEID field), SESET
    q     QSET/QSET1
    sa    CAEROi
    a     ASET/ASET1, CSUPEXT, superelement exterior DOFs
    c     CSET/CSET1, BNDREE/BNDFRE1
    b     BSET/BSET1, BNDFIX/BNDFIX1
    ap    ACCSPT
    rv    RVDOF/RVDOF1
    u1-u6 USET/USET1, SEUSET/SEUSET1

    * placed in the "m" set only if constraints are not specified in
      the basic coordinate system


    Set Name         Description
    --------         -----------
    s = sb + sg      all DOFs eliminated by single point constraints
    l = b + c + lm   the DOFs remaining after the reference DOFs are removed (DOF left over)
    t = l + r        the total set of physical boundary DOF for superelements
    a = t + q        the analysis set used in eigensolution
    d = a + e        the set used in dynamic analysis by the direct method
    f = a + o        unconstrained (free) structural DOFs
    fe = f + e       free DOFs plus extra DOFs
    n = f + s        all DOFs not constrained by multipoint constraints
    ne = n + e       all DOFs not constrained by multipoint constraints plus extra DOFs
    m = mp + mr      all DOFs eliminated by multipoint constraints
    g = n + m        all DOFs including scalar DOFs
    p = g + e        all physical DOFs including extra point DOFs
    ks = k + sa      the union of k and the re-used s-set (6 dof per grid)
    js = j + sa      the union of j and the re-used s-set (6 dof per grid)
    fr = o + l       statically independent set minus the statically determinate supports (fr = f – q – r)
    v = o + c + r    the set free to vibrate in dynamic reduction and component mode synthesis
    al = a – lm      a-set  without Lagrange multiplier DOFs
    dl = d – lm      d-set  without Lagrange multiplier DOFs
    gl = g – lm      g-set  without Lagrange multiplier DOFs
    ll = l – lm      l-set  without Lagrange multiplier DOFs
    nf = ne – lm     ne-set without Lagrange multiplier DOFs
    pl = p – lm      p-set  without Lagrange multiplier DOFs
    tl = t – lm      t-set  without Lagrange multiplier DOFs
    nl = n – lm      n-set  without Lagrange multiplier DOFs
    fl = f – lm      f-set  without Lagrange multiplier DOFs
    ff = fe – lm     fe-set without Lagrange multiplier DOFs


    Per MSC Nastran
    ---------------
    js = j + sa  the union of j and the re-used s-set (6 dof per grid)
    ks = k + sa  the union of k and the re-used s-set (6 dof per grid)

    [K]{x} = {F}
    a - active
    s - SPC
    [Kaa Kas]{xa} = {Fa}
    [Ksa Kss]{xs}   {Fs}
    """
    abs_kgg = np.abs(Kgg)
    col_kgg = abs_kgg.max(axis=0)
    row_kgg = abs_kgg.max(axis=1)
    # print(abs_kgg)
    # print(col_kgg)
    # print(row_kgg)
    # izero = np.where((col_kgg == 0.) & (row_kgg == 0))[0]
    if isinstance(Kgg, np.ndarray):
        ipositive = np.where((col_kgg > 0.0) | (row_kgg > 0))[0]
    elif issparse(Kgg):
        ipositive1 = col_kgg.todense().nonzero()[1]
        ipositive2 = row_kgg.todense().T.nonzero()[1]
        ipositive = np.union1d(ipositive1, ipositive2)
        if not isinstance(Kgg, csc_matrix):
            Kgg = Kgg.tocsc()
    else:
        raise NotImplementedError(type(Kgg))
    all_rows = np.arange(Kgg.shape[0], dtype=idtype)
    inegative = np.setdiff1d(all_rows, ipositive)
    apositive = aset[ipositive]
    sz_set = np.setdiff1d(aset, apositive)
    Kaa = Kgg[ipositive, :][:, ipositive]
    return Kaa, ipositive, inegative, sz_set


def guyan_reduction(matrix: np.ndarray, set1: np.ndarray, set2: np.ndarray):
    """
    https://en.wikipedia.org/wiki/Guyan_reduction

    [A11 A12]{d1} = {F1}
    [A21 A22]{d2} = {0}
    [A11]{d1} + [A12]{d2} = {F}
    [A21]{d1} + [A22]{d2} = {0}
    {d2} = -[A22]^-1[A21]{d1}
    [A11]{d1} + [A12](-[A22]^-1[A21]{d1}) = {F}
    ([A11] + ([A12]-[A22]^-1[A21]){d1} = {F}

    (AB)^-1 = B^-1 A^-1
    (cA)^-1 = A^-1/c
    (A)^-n = (A^-1)^n
    (A^-1)^T = (A^T)^-1

    https://math.stackexchange.com/questions/17776/inverse-of-the-sum-of-matrices
    (A + B)^-1 = A^-1 - 1/(1+g)A^-1 B A^-1; g=trace(B A^-1) where g != -1
               = A^-1 - A^-1*B*(A+B^-1)
    (A + B)^-1 = A^-1 + X
    (A^-1 + X)(A + B) = I
    A^-1 A + X A + A^-1 B + X B = I
    X(A + B) = A^-1 B
    X = -A^-1 B (A^-1 + B)^-1
    X = -A^-1 B (A^-1 + X)
    (I + A^-1 B) X = A^-1 B A^-1
    X = -(I + A^-1 B)^-1 (A^-1 B A^-1)
    """
    sets = [
        ["1", set1],
        ["2", set2],
    ]
    A = partition_matrix(matrix, sets)
    A11 = A["11"]
    A22 = A["22"]
    # A12 = A['12']
    A21 = A["21"]
    nA11 = A11.shape[0]

    # ([A11] + ([A12]-[A22]^-1[A21]){d1} = {F}
    A22m1 = np.linalg.inv(A22)
    # A12_A22m1_A21 = np.linalg.multi_dot([A11, A22m1, A21])
    T = np.vstack([np.eye(nA11), A22m1 @ A21])
    return np.linalg.multi_dot([T.T, A, T])
    # return A11 + A12_A22m1_A21


def _get_node_gridtype(model: BDF, idtype: str = "int32") -> NDArrayN2int:
    """
    Helper method for results post-processing

    Point type (per NX 10; OUG table; p.5-663):
    =1, GRID Point
    =2, Scalar Point
    =3, Extra Point
    =4, Modal
    =5, p-elements, 0-DOF
    -6, p-elements, number of DOF

    """
    node_gridtype = []
    # fdtype = 'float64'
    if model.grid.n:
        nodes = model.grid.node_id
        ones = np.ones(model.grid.n, dtype=idtype)
        node_gridtype.append(np.column_stack([nodes, ones]))

    if model.spoint.n:
        # raise NotImplementedError()
        spoint_id = model.spoint.spoint_id
        spoint_id.sort()
        twos = np.ones(model.spoint.n, dtype=idtype) * 2
        node_gridtype.append(np.column_stack([spoint_id, twos]))

    assert len(node_gridtype) > 0
    # print(node_gridtype)
    node_gridtype_array = np.vstack(node_gridtype)
    return node_gridtype_array


def get_aset(model: BDF) -> set[tuple[int, int]]:
    return model.aset.set_map


def get_bset(model: BDF) -> set[tuple[int, int]]:
    """creates the b-set"""
    return model.bset.set_map


def get_cset(model: BDF) -> set[tuple[int, int]]:
    """creates the c-set"""
    return model.cset.set_map


def get_omit_set(model: BDF) -> set[tuple[int, int]]:
    """creates the o-set"""
    return model.omit.set_map


def get_rset(model: BDF) -> set[tuple[int, int]]:
    """Creates the r-set from SUPORT/SUPORT1 cards.

    The r-set contains the reference DOFs used to define rigid-body motion
    in free-body (inertia relief) analysis.
    """
    rset_map: set[tuple[int, int]] = set()
    suport = model.suport
    if suport.n > 0:
        for nid, comp in zip(suport.node_id, suport.component):
            comp_str = str(comp)
            for compi in comp_str:
                if compi != '0':
                    rset_map.add((nid, int(compi)))
    return rset_map


def get_rset_bool(model: BDF,
                  dof_map: DOF_MAP,
                  ndof: int) -> np.ndarray:
    """Build a boolean r-set mask from SUPORT/SUPORT1 cards.

    Parameters
    ----------
    model : BDF
        The model.
    dof_map : DOF_MAP
        Maps (nid, dof) to index in g-set.
    ndof : int
        Total number of DOFs.

    Returns
    -------
    rset_b : (ndof,) bool array
        True at DOF indices belonging to the r-set (SUPORT reference DOFs).
    """
    rset_b = np.zeros(ndof, dtype='bool')
    rset_map = get_rset(model)
    for nid_dof in rset_map:
        if nid_dof in dof_map:
            rset_b[dof_map[nid_dof]] = True
    return rset_b


def partition_a_to_lr(
    Kaa: np.ndarray,
    Maa: np.ndarray,
    Fa: np.ndarray | None,
    aset_b: np.ndarray,
    rset_b: np.ndarray,
) -> dict[str, Any]:
    """Partition the a-set into l-set and r-set for SUPORT reduction.

    In Nastran's DOF partitioning: a = l + r, where r is defined by
    SUPORT/SUPORT1. The l-set (left-over) has rigid-body DOFs removed.

    Parameters
    ----------
    Kaa : (na, na) ndarray or sparse
        Stiffness in the analysis set.
    Maa : (na, na) ndarray
        Mass in the analysis set.
    Fa : (na,) ndarray or None
        Applied force in the analysis set. None for modes.
    aset_b : (ndof,) bool
        Boolean mask for the a-set in global DOFs.
    rset_b : (ndof,) bool
        Boolean mask for the r-set in global DOFs.

    Returns
    -------
    dict with keys:
        lset_local : integer index array into a-set for l-DOFs
        rset_local : integer index array into a-set for r-DOFs
        Kll, Klr, Krl, Krr : partitioned stiffness
        Mll, Mlr, Mrl, Mrr : partitioned mass
        Fl, Fr : partitioned force (None if Fa is None)
    """
    # Build local (within a-set) indices for l and r
    # aset indices in global DOFs
    a_indices = np.where(aset_b)[0]
    na = len(a_indices)

    # Which of the a-set DOFs are in the r-set?
    r_in_a = rset_b[a_indices]
    rset_local = np.where(r_in_a)[0]
    lset_local = np.where(~r_in_a)[0]

    # Convert Kaa to dense if sparse for partitioning
    if hasattr(Kaa, 'toarray'):
        Kaa_dense = Kaa.toarray()
    else:
        Kaa_dense = np.asarray(Kaa)

    K = partition_matrix(
        Kaa_dense, [("l", lset_local), ("r", rset_local)])
    M = partition_matrix(
        Maa, [("l", lset_local), ("r", rset_local)])

    result = {
        'lset_local': lset_local,
        'rset_local': rset_local,
        'Kll': K['ll'], 'Klr': K['lr'],
        'Krl': K['rl'], 'Krr': K['rr'],
        'Mll': M['ll'], 'Mlr': M['lr'],
        'Mrl': M['rl'], 'Mrr': M['rr'],
        'Fl': None, 'Fr': None,
    }
    if Fa is not None:
        result['Fl'] = Fa[lset_local]
        result['Fr'] = Fa[rset_local]
    return result


def get_qset(model: BDF) -> set[tuple[int, int]]:
    """creates the q-set"""
    return model.qset.set_map


def get_residual_structure(
    model: BDF, dof_map: DOF_MAP, fset: NDArrayNbool, idtype: str = "int32"
) -> NDArrayNbool:
    """gets the residual structure dofs"""
    asetmap = get_aset(model)
    bsetmap = get_bset(model)
    csetmap = get_cset(model)
    rsetmap = get_rset(model)
    qsetmap = get_qset(model)
    osetmap = get_omit_set(model)
    aset = apply_dof_map_to_set(asetmap, dof_map, idtype=idtype, use_ints=False)
    bset = apply_dof_map_to_set(bsetmap, dof_map, idtype=idtype, use_ints=False)
    cset = apply_dof_map_to_set(csetmap, dof_map, idtype=idtype, use_ints=False)
    rset = apply_dof_map_to_set(rsetmap, dof_map, idtype=idtype, use_ints=False)
    oset = apply_dof_map_to_set(osetmap, dof_map, idtype=idtype, use_ints=False)
    qset = apply_dof_map_to_set(qsetmap, dof_map, idtype=idtype, use_ints=False)

    # ndof = len(dof_map)
    # print('aset =', aset)
    # print('qset =', qset)
    # print('fset =', fset)
    # print('oset =', oset)

    aqo_set = np.vstack([aset, qset, oset])
    bcr_set = np.vstack([bset, cset, rset])
    # The a-set and o-set are created in the following ways:
    #    1. If only OMITi entries are present, then the o-set consists
    #       of degrees-of-freedom listed explicitly on OMITi entries.
    #       The remaining f-set degrees-of-freedom are placed in the
    #       b-set, which is a subset of the a-set.
    if np.any(oset) and ~np.any([aset, bset, cset, qset, rset]):
        # b = f - o
        bset = np.setdiff1d(fset, oset)
        aset = bset
        # 2. If ASETi or QSETi entries are present, then the a-set
        #    consists of all degrees-of-freedom listed on ASETi entries
        #    and any entries listing its subsets, such as QSETi, SUPORTi
        #    CSETi, and BSETi entries. Any OMITi entries are redundant.
        #    The remaining f-set degrees-of-freedom are placed in the
        #    o-set.
    elif np.any([aset, qset]):
        abcqr_set = np.vstack([aset, bset, cset, qset, rset])
        # print(abcqr_set)
        aset = np.sum(abcqr_set, axis=0).astype("bool")  # dim=2 -> row
        if np.any(oset):  # check no overlap with O set
            ao_set = aset | oset
            if np.any(ao_set):
                raise RuntimeError(
                    "OMITi entries cannot overlap with ASETi entries "
                    "or any ASET subsets, such as QSETi, "
                    "SUPORTi, CSETi, and BSETi entries."
                )

        oset = fset & ~aset  # assign remaining to O set
        # 3. If there are no ASETi, QSETi, or OMITi entries present but
        #    there are SUPORTi, BSETi, or CSETi entries present, then
        #    the entire f-set is placed in the a-set and the o-set is
        #    not created.
    elif not np.any(aqo_set) and np.any(bcr_set):
        aset = fset
        # 4. There must be at least one explicit ASETi, QSETi, or OMITi
        #    entry for the o-set to exist, even if the ASETi, QSETi, or
        #    OMITi entry is redundant. (related to item 3)
    else:
        # No model reduction - same as previous option
        aset = fset
    # Add ASET to residual structure TSET
    tset = aset & ~qset

    # Exclusive Degrees-of-freedom sets
    # ---------------------------------
    # m  # ([nGdof,1] logical) Degrees-of-freedom eliminated by multiple constraints
    # sb # ([nGdof,numSID] logical) Degrees-of-freedom eliminated by single-point constraints that are included in boundary conditions
    # sg # ([nGdof,1] logical) Degrees-of-freedom eliminated by single-point constraints that are specified on the PS field on node entries
    # o  # ([nGdof,1] logical) Degrees-of-freedom omitted by structural matrix partitioning
    # q  # ([nGdof,1] logical) Generalized degrees-of-freedom for dynamic reduction or component mode synthesis
    # r  # ([nGdof,1] logical) Reference degrees-of-freedom used to determine free body motion
    # c  # ([nGdof,1] logical) Degrees-of-freedom that are free during component mode synthesis or dynamic reduction
    # b  # ([nGdof,1] logical) Degrees-of-freedom fixed during component mode analysis or dynamic reduction
    # e  # ([nGdof,1] logical) extra degrees-of-freedom introduced in dynamic analysis
    # sa # Permanently constrained aerodynamic degrees-of-freedom
    # k  # Aerodynamic degrees-of-freedom

    # Nonexclusive Degrees-of-freedom sets
    # ------------------------------------
    # s  # ([nGdof,numSID] logical) [sb + sg] Degrees-of-freedom eliminated by single point constraints
    # l  # ([nGdof,1] logical) [b + c] Structural degrees-of-freedom remaining after the reference degrees-of-freedom are removed (degrees-of-freedom left over)
    #    # ([nGdof,1] logical) [l + r] Total set of physical boundary degrees-of-freedom for superelements
    #    # ([nGdof,1] logical) [t + q] Set assembled in superelement analysis
    # d  # ([nGdof,1] logical) [a + e] Set used in dynamic analysis by the direct method
    #    # ([nGdof,1] logical) [a + o] Unconstrained (free) structural degrees-of-freedom
    # fe # ([nGdof,1] logical) [f + e] Free structural degrees-of-freedom plus extra degrees-of-freedom
    #    # ([nGdof,1] logical) [f + s] Degrees-of-freedom not constrained by multipoint constraints
    # ne % ([nGdof,1] logical) [n + e] Structural degrees-of-freedom not constrained by multipoint constraints plus extra degrees-of-freedom
    # g = true(nGdof,1) [n + m] All structural degrees-of-freedom including scalar degrees-of-freedom
    # p = [g + e] Physical degrees-of-freedom
    # ps = [p + sa] Physical and constrained (SPCi) aerodynamic degrees-of-freedom
    # pa = [ps + k] Physical set for aerodynamics
    # fr = [f ? q ? r] Statically independent set minus the statically determinate supports
    # v = [o + c + r] Set free to vibrate in dynamic reduction and component mode synthesis
    return aset, tset


def apply_dof_map_to_set(
    set_map, dof_map: DOF_MAP, idtype: str = "int32", use_ints: bool = True
) -> NDArrayNbool:
    """changes a set defined in terms of (nid, comp) into an array of integers"""
    if use_ints:
        ndof = len(set_map)
        aset = np.full(ndof, 0, dtype=idtype)
        for i, dofi in enumerate(set_map):
            aset[i] = dof_map[dofi]
        aset.sort()
    else:
        ndof = len(dof_map)
        aset = np.full(ndof, 0, dtype="bool")
        for dofi in set_map:
            i = dof_map[dofi]
            aset[i] = True
    return aset


def xg_to_xb(
    model, xg: NDArrayNfloat, ngrid: int, ndof_per_grid: int, inplace: bool = True
) -> NDArrayNfloat:
    assert isinstance(xg, np.ndarray)
    str(ngrid)

    xb = xg
    if not inplace:
        xb = copy.deepcopy(xb)

    nids = model._type_to_id_map["GRID"]
    for i, nid in enumerate(nids):
        node = model.nodes[nid]
        if node.cd:
            model.log.debug(f"node {nid} has a CD={node.cd}")
            cd_ref = node.cd_ref
            T = cd_ref.beta_n(n=2)
            i1 = i * ndof_per_grid
            i2 = (i + 1) * ndof_per_grid
            xi = xg[i1:i2]
            xb[i1:i2] = xi @ T  # TODO: verify the transform; I think it's right
    return xb


def write_oload(
    Fb: NDArrayNfloat,
    dof_map: DOF_MAP,
    isubcase: int,
    ngrid: int,
    ndof_per_grid: int,
    f06_file: TextIO,
    page_stamp: str,
    page_num: int,
    log: SimpleLogger,
) -> int:
    """writes the OLOAD RESULTANT table"""
    str(ngrid)
    str(dof_map)
    fxyz_mxyz = Fb[: ngrid * ndof_per_grid].reshape(ngrid, ndof_per_grid)
    fxyz_mxyz_sum = fxyz_mxyz.sum(axis=0)

    f06_file.write(
        " *** USER INFORMATION MESSAGE 7310 (VECPRN)\n"
        "     ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM WILL BE USED AS REFERENCE LOCATION.\n"
        "     RESULTANTS ABOUT ORIGIN OF SUPERELEMENT BASIC COORDINATE SYSTEM IN "
        "SUPERELEMENT BASIC SYSTEM COORDINATES.\n"
    )
    oload = Resultant("OLOAD", fxyz_mxyz_sum, isubcase)
    log.info(f"OLOAD {fxyz_mxyz_sum}")
    page_num = oload.write_f06(f06_file, page_stamp, page_num)
    return page_num + 1


def _write_spcforce_resultant(
    f06_file: TextIO,
    fspc: NDArrayNfloat,
    ngrid: int,
    ndof_per_grid: int,
    isubcase: int,
    page_stamp: str,
    page_num: int,
    log: SimpleLogger,
) -> None:
    """Write SPCFORCE RESULTANT table to F06."""
    fxyz_mxyz = fspc[:ngrid * ndof_per_grid].reshape(ngrid, ndof_per_grid)
    fxyz_mxyz_sum = fxyz_mxyz.sum(axis=0)
    spc_resultant = Resultant("SPCFORCE", fxyz_mxyz_sum, isubcase)
    log.info(f"SPCFORCE RESULTANT {fxyz_mxyz_sum}")
    spc_resultant.write_f06(f06_file, page_stamp, page_num)


def solve(
    Kaa: lil_matrix,
    Fa_solve: np.ndarray,
    aset: np.ndarray,
    log: SimpleLogger,
    idtype: str = "int32",
):
    """solves [K]{u} = {F}"""
    log.info("starting solve")
    Kaa_, ipositive, inegative, unused_sz_set = remove_rows(Kaa, aset, idtype=idtype)
    if len(ipositive) == 0:
        log.error("empty Kaa")
    # if np.linalg.det(Kaa) == 0.:
    # log.error('singular Kaa')

    isolve = len(ipositive)
    if isolve == 0:
        log.error(f"  ipositive = {ipositive}")
        log.error(f"  Kaa_ = {Kaa_}")
        # print(Kaa_)
        raise RuntimeError("no residual structure found")

    # Maa_ = Maa[ipositive, :][:, ipositive]
    # print(f'Fg = {Fg}')
    # print(f'Fa = {Fa}')
    # print(f'Fs = {Fs}')

    Fa_ = Fa_solve[ipositive]
    # [A]{x} = {b}
    # [Kaa]{x} = {F}
    # {x} = [Kaa][F]
    # print(f'Kaa:\n{Kaa}')
    # print(f'Fa: {Fa}')

    Kaa_dense = Kaa_.toarray()
    np.set_printoptions(precision=0, linewidth=100)
    log.debug(f"  Kaas_:\n{Kaa_dense}")
    log.debug(f"  Kaa_:\n{Kaa_}")
    log.debug(f"  Fa_: {Fa_}")

    backend = get_solver()
    log.debug(f"  solver backend: {backend.name}")
    xa_ = np.linalg.solve(Kaa_dense, Fa_)
    xas_ = backend.spsolve(Kaa_, Fa_)
    sparse_error = np.linalg.norm(xa_ - xas_)
    if sparse_error > 1e-12:
        log.warning(f"  sparse_error = {sparse_error}")
    log.info("finished solve")
    return xas_, ipositive, inegative


def build_Mbb(
    model: BDF, subcase: Subcase, dof_map: DOF_MAP, ndof: int, fdtype="float64"
) -> NDArrayNNfloat:
    """builds the mass matrix in the basic frame, [Mbb]"""
    log = model.log
    log.info("starting build_Mbb")
    wtmass = 1.0
    from pyNastran.dev.bdf_vectorized3.solver.build_stiffness import _COOAccumulator
    _coo_m = _COOAccumulator(ndof)
    str(model)
    str(subcase)
    # no_mass = {
    #'CELAS1', 'CELAS2', 'CELAS3', 'CELAS4',
    #'CDAMP1', 'CDAMP2', 'CDAMP3', 'CDAMP4',
    # }

    mass_rod_2x2 = (
        np.array(
            [
                [2, 0, 1, 0],
                [0, 2, 0, 1],
                [1, 0, 2, 0],
                [0, 1, 0, 2],
            ],
            dtype="float64",
        )
        / 3.0
    )
    mass_tri = (
        np.array(
            [
                [4, 0, 2, 0, 1, 0],
                [0, 4, 0, 2, 0, 1],
                [2, 0, 4, 0, 2, 0],
                [0, 2, 0, 4, 0, 2],
                [1, 0, 2, 0, 4, 0],
                [0, 1, 0, 2, 0, 4],
            ],
            dtype="float64",
        )
        / 6.0
    )
    # mass_quad_1x1 = np.array([
    # [1, 0, 1, 0, 1, 0, 1, 0],
    # [0, 1, 0, 1, 0, 1, 0, 1],
    # [1, 0, 1, 0, 1, 0, 1, 0],
    # [0, 1, 0, 1, 0, 1, 0, 1],
    # [1, 0, 1, 0, 1, 0, 1, 0],
    # [0, 1, 0, 1, 0, 1, 0, 1],
    # [1, 0, 1, 0, 1, 0, 1, 0],
    # [0, 1, 0, 1, 0, 1, 0, 1],
    # ], dtype='float64') / 36.

    mass_quad_2x2 = (
        np.array(
            [
                [4, 0, 2, 0, 1, 0, 2, 0],
                [0, 4, 0, 2, 0, 1, 0, 2],
                [2, 0, 4, 0, 2, 0, 1, 0],
                [0, 2, 0, 4, 0, 2, 0, 1],
                [1, 0, 2, 0, 4, 0, 2, 0],
                [0, 1, 0, 2, 0, 4, 0, 2],
                [2, 0, 1, 0, 2, 0, 4, 0],
                [0, 2, 0, 1, 0, 2, 0, 4],
            ],
            dtype="float64",
        )
        / 36.0
    )

    mass_total = 0.0
    if model.conm1:
        nid = elem.nid
        nid_ref = elem.nid_ref
        if nid_ref.type == "GRID":
            i1 = dof_map[(nid, 1)]

            # TODO: support CID
            if nid_ref.cd != elem.cid:
                log.warning(
                    f"  CONM1 eid={eid} nid={nid} CD={nid_ref.cd} to cid={elem.cid} is not supported"
                )
        else:  # pragma: no cover
            print(elem.get_stats())
            raise NotImplementedError(elem)
        _coo_m.add_matrix(list(range(i1, i1 + 6)), elem.mass_matrix)

    # PARAM,COUPMASS: 1 => consistent mass, -1 or absent => lumped
    use_consistent_shells = False
    if hasattr(model, "params") and "COUPMASS" in model.params:
        coupmass_val = model.params["COUPMASS"].values[0]
        if coupmass_val >= 1:
            use_consistent_shells = True

    # Consistent mass matrices for shells (m/36 * [4,2,1,2;...])
    mass_quad_consistent = (
        np.array(
            [
                [4, 2, 1, 2],
                [2, 4, 2, 1],
                [1, 2, 4, 2],
                [2, 1, 2, 4],
            ],
            dtype="float64",
        )
        / 36.0
    )
    mass_tri_consistent = (
        np.array(
            [
                [2, 1, 1],
                [1, 2, 1],
                [1, 1, 2],
            ],
            dtype="float64",
        )
        / 12.0
    )

    # has possibility of mass
    has_mass = False

    for elem in model.element_cards:
        etype = elem.type
        if etype in NO_MASS:
            continue
        if elem.n == 0:
            continue

        has_mass = True
        if etype == "CONM2":
            mass_total = conm2_fill_Mbb(model, mass_total, _coo_m, dof_map)

        elif etype in ["CROD", "CONROD", "CTUBE"]:
            # verified
            mass = elem.mass()
            if mass == 0.0:
                log.warning(f"  no mass for {etype} eid={eid}")
                continue

            nids1 = elem.nodes[:, 0]
            nids2 = elem.nodes[:, 1]
            for nid1, nid2 in zip(nids1, nids2):
                i1 = dof_map[(nid1, 1)]
                j1 = dof_map[(nid2, 1)]
                ii = [i1, i1 + 1, j1, j1 + 1]
                _coo_m.add_matrix(ii, mass_rod_2x2 * mass)
        elif etype in {"CBAR", "CBEAM"}:
            from pyNastran.dev.bdf_vectorized3.solver.beam import (
                consistent_mass,
                lumped_mass,
                beam_transform,
            )
            # PARAM,COUPMASS,1 => consistent; default (0 or absent) => lumped
            use_consistent = False
            if hasattr(model, 'params') and 'COUPMASS' in model.params:
                coupmass_val = model.params['COUPMASS'].values[0]
                if coupmass_val >= 1:
                    use_consistent = True

            area = elem.area()
            inertia = elem.inertia()
            xyz1, xyz2 = elem.get_xyz()
            lengths = np.linalg.norm(xyz2 - xyz1, axis=1)
            v, ihat_arr, yhat_arr, zhat_arr, wa_arr, wb_arr = elem.get_axes(xyz1, xyz2)
            k_arr = elem.k()
            e_g_nus = elem.e_g_nu()

            elem_masses = elem.mass()
            mass_total += elem_masses.sum()
            mass_per_length_total = elem_masses / lengths

            for (nid1, nid2), areai, inertiai, lengthi, ki, e_g_nu, ihati, jhati, khati, mpl in zip(
                elem.nodes, area, inertia, lengths, k_arr, e_g_nus,
                ihat_arr, yhat_arr, zhat_arr, mass_per_length_total,
            ):
                i1_inertia, i2_inertia, i12, j = inertiai
                e, g, nu = e_g_nu
                k1, k2 = ki

                rho_eff = mpl / areai if areai > 0 else 0.0

                if use_consistent:
                    Me = consistent_mass(
                        areai, lengthi, rho_eff,
                        i1_inertia, i2_inertia, j,
                        k1, k2, nsm=0.0,
                    )
                else:
                    Me = lumped_mass(
                        areai, lengthi, rho_eff,
                        i1_inertia, i2_inertia, j,
                        nsm=0.0,
                    )

                Teb = beam_transform(ihati, jhati, khati)
                M_basic = Teb.T @ Me @ Teb

                gi1 = dof_map[(nid1, 1)]
                gi2 = dof_map[(nid2, 1)]
                n_ijv = [
                    gi1, gi1 + 1, gi1 + 2, gi1 + 3, gi1 + 4, gi1 + 5,
                    gi2, gi2 + 1, gi2 + 2, gi2 + 3, gi2 + 4, gi2 + 5,
                ]
                _coo_m.add_matrix(n_ijv, M_basic)
        elif etype == "CTRIA3":
            masses = elem.mass()
            if masses.sum() == 0.0:
                log.warning(f"  no mass for CTRIA3 eid={elem.element_id}")
                continue
            for (nid1, nid2, nid3), massi in zip(elem.nodes, masses):
                if massi == 0.0:
                    continue
                if use_consistent_shells:
                    # Consistent: m * [2,1,1;1,2,1;1,1,2]/12
                    nids_e = [nid1, nid2, nid3]
                    M_consist = mass_tri_consistent * massi
                    for dof_offset in range(3):
                        ii = [dof_map[(n, 1)] + dof_offset for n in nids_e]
                        _coo_m.add_matrix(ii, M_consist)
                else:
                    # Lumped: m/3 per node on Tx, Ty, Tz
                    m_node = massi / 3.0
                    for nid in [nid1, nid2, nid3]:
                        i1 = dof_map[(nid, 1)]
                        _coo_m.add_scalar(i1, i1, m_node)
                        _coo_m.add_scalar(i1 + 1, i1 + 1, m_node)
                        _coo_m.add_scalar(i1 + 2, i1 + 2, m_node)
            # Mbb[i1, i1] = Mbb[i1+1, i1+1] = Mbb[i1+2, i1+2] = \
            # Mbb[i2, i2] = Mbb[i2+1, i2+1] = Mbb[i2+2, i2+2] = \
            # Mbb[i3, i3] = Mbb[i3+1, i3+1] = Mbb[i3+2, i3+2] = mass / 3
        elif etype == "CQUAD4":
            masses = elem.mass()
            if masses.sum() == 0.0:
                log.warning(f"  no mass for CQUAD4 eid={elem.element_id}")
                continue
            for (nid1, nid2, nid3, nid4), massi in zip(elem.nodes, masses):
                if massi == 0.0:
                    continue
                if use_consistent_shells:
                    # Consistent: m * [4,2,1,2;2,4,2,1;1,2,4,2;2,1,2,4]/36
                    # Applied to each translational DOF (Tx, Ty, Tz) independently
                    nids_e = [nid1, nid2, nid3, nid4]
                    M_consist = mass_quad_consistent * massi
                    for dof_offset in range(3):
                        ii = [dof_map[(n, 1)] + dof_offset for n in nids_e]
                        _coo_m.add_matrix(ii, M_consist)
                else:
                    # Lumped: m/4 per node on Tx, Ty, Tz
                    m_node = massi / 4.0
                    for nid in [nid1, nid2, nid3, nid4]:
                        i1 = dof_map[(nid, 1)]
                        _coo_m.add_scalar(i1, i1, m_node)
                        _coo_m.add_scalar(i1 + 1, i1 + 1, m_node)
                        _coo_m.add_scalar(i1 + 2, i1 + 2, m_node)
            # if 0:  # pragma: no cover
            # mass4 = mass / 9. # 4/36
            # mass2 = mass / 18. # 2/36
            # mass1 = mass / 36.
            # print(mass1, mass2, mass4)
            # Mbb[i1, i1] += mass4
            # Mbb[i1+1, i1+1] += mass4
            # Mbb[i2, i2] += mass4
            # Mbb[i2+1, i2+1] += mass4
            # Mbb[i3, i3] += mass4
            # Mbb[i3+1, i3+1] += mass4
            # Mbb[i4, i4] += mass4
            # Mbb[i4+1, i4+1] += mass4

            # Mbb[i1, i3] += mass1
            # Mbb[i1, i2] += mass2
            # Mbb[i1, i4] += mass2
            # Mbb[i3, i1] += mass1
            # Mbb[i2, i1] += mass2
            # Mbb[i4, i1] += mass2

            # Mbb[i1+1, i3+1] += mass1
            # Mbb[i1+1, i2+1] += mass1
            # Mbb[i1+1, i4+1] += mass2
            # Mbb[i3+1, i1+1] += mass1
            # Mbb[i2+1, i1+1] += mass2
            # Mbb[i4+1, i1+1] += mass2

            # Mbb[i2, i4] += mass1
            # Mbb[i2+1, i4+1] += mass1
            # Mbb[i2, i3] += mass2
            # Mbb[i2+1, i3+1] += mass2

            # Mbb[i4, i2] += mass1
            # Mbb[i4+1, i2+1] += mass1
            # Mbb[i3, i2] += mass2
            # Mbb[i3+1, i2+1] += mass2

            # Mbb[i3, i4] += mass2
            # Mbb[i3+1, i4+1] = mass2
            # Mbb[i4, i3] += mass2
            # Mbb[i4+1, i3+1] += mass2

            # Mbb[i2+1, i3+1] = 1
            # Mbb[i2+1, i2+1] = Mbb[i1+1, i4+1] = 2

            # Mbb[i2, i1+1] = Mbb[i3, i1+1] = Mbb[i4, i1+1] = 1
            # Mbb[i1, i2+1] = Mbb[i3, i2+1] = Mbb[i4, i2+1] = 1
            # Mbb[i1, i3+1] = Mbb[i2, i3+1] = Mbb[i4, i3+1] = 1
            # Mbb[i1, i4+1] = Mbb[i2, i4+1] = Mbb[i3, i4+1] = 1
            # print(Mbb)
            # print(Mbb[ii, :][:, ii])
        elif etype == "CSHEAR":
            masses = elem.mass()
            if masses.sum() == 0.0:
                continue
            for (nid1, nid2, nid3, nid4), massi in zip(elem.nodes, masses):
                if massi == 0.0:
                    continue
                i1 = dof_map[(nid1, 1)]
                i2 = dof_map[(nid2, 1)]
                i3 = dof_map[(nid3, 1)]
                i4 = dof_map[(nid4, 1)]
                ii = [
                    i1,
                    i1 + 1,
                    i2,
                    i2 + 1,
                    i3,
                    i3 + 1,
                    i4,
                    i4 + 1,
                ]
                _coo_m.add_matrix(ii, mass_quad_2x2 * massi)
        elif etype in {"CHEXA", "CTETRA", "CPENTA"}:
            pass  # handled below
        else:  # pragma: no cover
            print(elem.get_stats())
            raise NotImplementedError(elem)

    # Solid elements (CHEXA, CTETRA, CPENTA) — consistent mass
    from pyNastran.dev.bdf_vectorized3.solver.solids import build_mbb_solids
    mass_total += build_mbb_solids(model, _coo_m, dof_map)

    # Convert COO accumulator to sparse CSC
    Mbb = _coo_m.to_csc()

    if wtmass != 1.0:
        Mbb *= wtmass

    has_special_points = "SPOINT" in model.card_count or "EPOINT" in model.card_count
    is_all_grids = not has_special_points
    unused_can_dof_slice = is_all_grids and not has_mass
    if Mbb.nnz > 0:
        i = np.arange(0, ndof).reshape(ndof // 6, 6)[:, :3].ravel()
        massi = sum(Mbb[ii, ii] for ii in i)
        log.info(f"finished build_Mbb; M={massi:.6g}; mass_total={mass_total:.6g}")
    else:
        return None
    return Mbb


def grid_point_weight(model: BDF, Mbb, dof_map: DOF_MAP, ndof: int):
    str(dof_map)
    str(ndof)
    z = np.zeros((3, 3), dtype="float64")
    nnodes = model.grid.n
    nspoints = model.spoint.n
    # nepoints = len(model.epoints)
    ndof_grid = 6 * nnodes
    ndof2 = 6 * nnodes + nspoints
    # print(f'ndof={ndof} ndof2={ndof2}')
    D = np.zeros((ndof2, 6), dtype="float64")
    D[ndof_grid:, 0] = 1
    # print('D.shape ', D.shape)
    inid = 0
    # print(model.nodes)
    # print(f'D.shape = {D.shape}')
    reference_point = 0  #  model.get_param('GRDPNT', 0)
    if reference_point == 0:
        dxyz = np.zeros(3, dtype=Mbb.dtype)
    else:
        dxyz = model.nodes[reference_point].get_position()

    coord = model.coord
    xyz_cid0 = model.grid.xyz_cid0()
    cds = model.grid.cd
    cids = coord.coord_id
    # icds = coord.index(cds, assume_sorted=False, inverse=False)
    beta = coord.j
    for xyz0, cd in zip(xyz_cid0, cds):
        # coord = model.coord.s
        # print(f'nid={nid}')

        # TODO: what about otuput coordinate frames; specifically cylindrical and spherical frames?
        xi, yi, zi = xyz0 - dxyz
        Tr = np.array(
            [
                [0, zi, -yi],
                [-zi, 0, xi],
                [yi, -xi, 0],
            ],
            dtype="float64",
        )
        # cd_ref = node.cd_ref
        beta = coord.xyz_to_global_transform[cd]
        Ti = beta
        TiT_Tr = Ti.T @ Tr
        d = np.block(
            [
                [Ti.T, TiT_Tr],
                [z, Ti.T],
            ]
        )
        # print(f'd.shape={d.T.shape}')
        # print(f'd={d}')
        # print(f'D[:, {inid*6}:{(inid+1)*6}]')
        D[inid * 6 : (inid + 1) * 6, :] = d.T
        inid += 1

    # Location Basic System (CP Fields) Mass Global System (CD Fields)
    # Grid ID xb yb zb xCD yCD zCD
    # 1       0   0 0  2   3   5
    # 2       1   0 0  2   3   5
    # 3       0.5 1 0  2   3   5
    # 4       0.5 0 1  2   3   5
    # print(f'Dt.shape = {D.T.shape}')
    # print(f'Mbb.shape = {Mbb.shape}')
    # print(f'D.shape = {D.shape}')
    # print(f'D.T =\n{D.T}')
    M0 = D.T @ Mbb @ D
    return reference_point, M0


def dof_map_to_tr_set(dof_map, ndof: int) -> tuple[NDArrayNbool, NDArrayNbool]:
    """creates the translation/rotation sets"""
    trans_set = np.zeros(ndof, dtype="bool")
    rot_set = np.zeros(ndof, dtype="bool")
    for (unused_nid, dof), idof in dof_map.items():
        if dof in [0, 1, 2, 3]:
            trans_set[idof] = True
        else:
            rot_set[idof] = True
    return trans_set, rot_set


def ps_to_sg_set(ndof: int, ps: list[int]):
    """creates the SPC on the GRID (PS-field) set, {sg}"""
    # all DOFs are initially assumed to be active
    sg_set = np.ones(ndof, dtype="bool")

    # False means it's constrained
    sg_set[ps] = False
    return sg_set


def _has_rigid_elements(model: BDF) -> bool:
    """Check if the model has any rigid elements (RBE2, RBE3, etc.)."""
    for attr in ("rbe2", "rbe3", "rbar", "rbar1", "rbe1", "rrod"):
        elem = getattr(model, attr, None)
        if elem is not None and elem.n > 0:
            return True
    return False


def autospc_n_set(
    Knn,
    ndof_n: int,
    log,
    eps_ratio: float = 1.0e-8,
) -> np.ndarray:
    """Detect and remove singular DOFs from the N-set stiffness (AUTOSPC).

    After MPC/rigid element reduction, some DOFs in the N-set may be
    unconstrained (zero or near-zero stiffness). This function identifies
    them so they can be added to the SPC set.

    Parameters
    ----------
    Knn : sparse or dense matrix (ndof_n, ndof_n)
        N-set stiffness matrix.
    ndof_n : int
        Size of the N-set.
    log : logger
        Logging object.
    eps_ratio : float
        A diagonal entry is flagged as singular if
        |Kii| / max(|diag(K)|) < eps_ratio.

    Returns
    -------
    autospc_dofs : (n_singular,) int array
        Indices of singular DOFs in the N-set (to be constrained).
    """
    if issparse(Knn):
        diag = np.abs(np.array(Knn.diagonal()).ravel())
    else:
        diag = np.abs(np.diag(Knn))

    max_diag = diag.max() if len(diag) > 0 else 1.0
    if max_diag == 0.0:
        return np.arange(ndof_n, dtype="int32")

    singular = np.where(diag / max_diag < eps_ratio)[0].astype("int32")

    if len(singular) > 0:
        log.warning(
            f"  AUTOSPC: {len(singular)} near-singular DOFs detected in N-set "
            f"(ratio < {eps_ratio:.1e}), constraining them"
        )
    return singular


def _get_dof_map(model: BDF) -> tuple[DOF_MAP, list[int]]:
    """helper method for ``get_static_force_vector_from_subcase_id``"""
    i = 0
    dof_map = {}
    spoints: list[int] = []
    ps: list[int] = []
    for nid, psi in zip(model.grid.node_id, model.grid.ps):
        for dof in range(1, 7):
            dof_map[(nid, dof)] = i
            i += 1
        if psi != 0:
            asf
            for psii in str(psi):
                nid_dof = (nid, int(psii))
                j = dof_map[nid_dof]
                ps.append(j)

    # we want the GRID points to be first
    assert len(spoints) == 0, spoints

    if model.spoint.n:
        for nid in model.spoint.spoint_id:
            key = (nid, 0)
            if key not in dof_map:
                dof_map[key] = i
                i += 1
    assert len(dof_map) > 0
    return dof_map, ps


def get_ndof(model: BDF, subcase: Subcase) -> tuple[int, int, int]:
    ndof_per_grid = 6
    if "HEAT" in subcase:
        ndof_per_grid = 1
    ngrid = model.grid.n
    nspoint = model.spoint.n
    # nepoint = model.epoint.n
    nepoint = 0
    ndof = ngrid * ndof_per_grid + nspoint + nepoint
    # print(f'ngrid={ngrid} nspoint={nspoint}')
    assert ndof > 0, model.card_count
    return ngrid, ndof_per_grid, ndof


def conm2_fill_Mbb(
    model: BDF, mass_total: float, _coo_m, dof_map: dict[tuple[int, int], int]
) -> float:
    eye3 = np.eye(3, dtype="float64")
    conm2 = model.conm2
    log = model.log
    # mass = elem.Mass()
    # nid = elem.nid
    # nid_ref = elem.nid_ref
    inid = model.grid.index(conm2.node_id)
    cds = model.grid.cd[inid]
    for eid, nid, cd, cid, mass, elem_x, elem_i in zip(
        conm2.element_id,
        conm2.node_id,
        cds,
        conm2.coord_id,
        conm2.mass(),
        conm2.xyz_offset,
        conm2.inertia,
    ):
        i1 = dof_map[(nid, 1)]
        if cd != cid:
            log.warning(f"  CONM2 eid={eid} nid={nid} CD={cd} to CONM2 cid={cid} is not supported")
        # Mbb[i1, i1] = mass
        # Mbb[i1+1, i1+1] = mass
        # Mbb[i1+2, i1+2] = mass
        # TODO: support CID
        I11, I21, I22, I31, I32, I33 = elem_i
        x1, x2, x3 = elem_x
        mxx = (
            np.array(
                [
                    [x1 * x1, -x1 * x2, -x1 * x3],
                    [-x2 * x1, x2 * x2, -x2 * x3],
                    [-x3 * x1, x3 * x2, x3 * x3],
                ]
            )
            * mass
        )
        Tr = np.array(
            [
                [0, x3, -x2],
                [-x3, 0, x1],
                [x2, -x1, 0],
            ],
            dtype="float64",
        )
        mx = Tr * mass
        I = (
            np.array(
                [
                    [I11, -I21, I31],
                    [-I21, I22, -I32],
                    [-I31, -I32, I33],
                ]
            )
            + mxx
        )

        # [mass, 01, 02, 03, mass * X3, -mass * X2]
        # [10, mass, 12, -mass * X3, 14, mass * X1]
        # [20, 21, mass, mass * X2, -mass * X1, 25]
        # [30, -mass * X3, mass * X2,        I11 + mass * X2 * X2 + mass * X3 * X3, -I21 - mass * X2 * X1,                  -I31 - mass * X3 * X1]
        # [mass * X3, 41, -mass * X1,       -I21 - mass * X2 * X1,                   I22 + mass * X1 * X1 + mass * X3 * X3, -I32 - mass * X3 * X2]
        # [-mass * X2, mass * X1, 52,       -I31 - mass * X3 * X1,                  -I32 - mass * X3 * X2,                   I33 + mass * X2 * X2 + mass * X1 * X1]

        M6 = np.zeros((6, 6), dtype='float64')
        M6[:3, :3] = eye3 * mass
        M6[:3, 3:] = mx
        M6[3:, :3] = mx.T
        M6[3:, 3:] = I
        _coo_m.add_matrix(list(range(i1, i1 + 6)), M6)
        mass_total += mass
    return mass_total


def _run_modes(
    model: BDF,
    subcase: Subcase,
    ngrid: int,
    ndof_per_grid: int,
    ndof: int,
    node_gridtype: np.ndarray,
    dof_map: dict,
    Kgg_override=None,
    GMN=None,
    mset: np.ndarray | None = None,
    idtype: str = "int32",
    fdtype: str = "float64",
):
    log = model.log
    out = {}
    if Kgg_override is not None:
        if issparse(Kgg_override):
            Kgg = Kgg_override
        else:
            Kgg = csc_matrix(Kgg_override[:ndof, :ndof])
    else:
        Kgg = build_Kgg(model, dof_map, ndof, ngrid, ndof_per_grid, idtype="int32", fdtype=fdtype)
    out["Kgg"] = Kgg

    Mbb = build_Mbb(model, subcase, dof_map, ndof, fdtype=fdtype)

    Mgg = Kbb_to_Kgg(model, Mbb, ngrid, ndof_per_grid)
    out["Mgg"] = Mgg
    del Mbb

    # ----------------------MPC reduction (g -> n)----------------------
    if GMN is not None:
        Knn = GMN.T @ Kgg @ GMN
        Mnn = GMN.T @ Mgg @ GMN
        # Work in n-set
        gset_all = np.arange(ndof, dtype="int32")
        nset_indices = np.setdiff1d(gset_all, mset)
        n_ndof = len(nset_indices)

        # Map g-set SPC to n-set
        sset_g, sset_b_g, xg = _build_xg(model, dof_map, ndof, subcase)
        g_to_n = np.full(ndof, -1, dtype="int32")
        g_to_n[nset_indices] = np.arange(n_ndof, dtype="int32")

        sset_b = np.zeros(n_ndof, dtype="bool")
        for s_idx in sset_g:
            n_idx = g_to_n[s_idx]
            if n_idx >= 0:
                sset_b[n_idx] = True
        sset = sset_b

        gset = np.ones(n_ndof, dtype="bool")
        aset = gset & ~sset_b
        out["aset"] = aset
        out["xg"] = xg

        xn = np.zeros(n_ndof, dtype=fdtype)
        xa, xs = partition_vector2(xn, [["a", aset], ["s", sset_b]])

        M = partition_matrix(Mnn, [["a", aset], ["s", sset_b]])
        Maa = M["aa"]
        out["Maa"] = Maa

        K = partition_matrix(Knn, [["a", aset], ["s", sset_b]])
        Kaa = K["aa"]
        out["Kaa"] = Kaa

        ndof_work = n_ndof
        Kgg_work = Knn
        Mgg_work = Mnn
    else:
        gset = np.ones(ndof, dtype="bool")
        sset, sset_b, xg = _build_xg(model, dof_map, ndof, subcase)
        aset = gset & ~sset_b  # a = g-s
        out["aset"] = aset

        xa, xs = partition_vector2(xg, [["a", aset], ["s", sset_b]])
        out["xg"] = xg

        M = partition_matrix(Mgg, [["a", aset], ["s", sset_b]])
        Maa = M["aa"]
        out["Maa"] = Maa

        K = partition_matrix(Kgg, [["a", aset], ["s", sset_b]])
        Kaa = K["aa"]
        out["Kaa"] = Kaa

        ndof_work = ndof
        Kgg_work = Kgg
        Mgg_work = Mgg
    # Kss = K['ss']
    # Kas = K['as']
    # Ksa = K['sa']

    # [Kaa]{xa} + [Kas]{xs} = {Fa}
    # [Ksa]{xa} + [Kss]{xs} = {Fs}

    # {xa} = [Kaa]^-1 * ({Fa} - [Kas]{xs})
    # {Fs} = [Ksa]{xa} + [Kss]{xs}

    # --- R-set partitioning (SUPORT reduction) ---
    # If SUPORT cards exist, partition a = l + r.
    # Solve eigenvalue problem on the l-set (flexible DOFs).
    rset_b = get_rset_bool(model, dof_map, ndof)
    has_rset = np.any(rset_b) and GMN is None

    neigenvalue, norm_str = get_real_eigenvalue_method(model, subcase)

    if has_rset:
        # r-set within the a-set
        a_indices = np.where(aset)[0]
        r_in_a = rset_b[a_indices]
        lset_local = np.where(~r_in_a)[0]
        rset_local = np.where(r_in_a)[0]
        nr = len(rset_local)
        nl = len(lset_local)
        log.info(f"SUPORT r-set: {nr} DOFs, l-set: {nl} DOFs")

        if hasattr(Kaa, 'toarray'):
            Kaa_dense = Kaa.toarray()
        else:
            Kaa_dense = np.asarray(Kaa)

        Kll = Kaa_dense[np.ix_(lset_local, lset_local)]
        Mll = Maa[np.ix_(lset_local, lset_local)]

        # Remove zero rows/columns (AUTOSPC equivalent for l-set)
        Kll_sp = csc_matrix(Kll)
        Kll_, ipositive, unused_ineg, unused_sz = remove_rows(
            Kll_sp, np.ones(nl, dtype='bool'))
        Mll_ = Mll[ipositive, :][:, ipositive]
        ndof_ = Kll_.shape[0]

        if issparse(Mll_):
            Mll_ = Mll_.toarray()
        eigenvalue, xl_ = solve_eigenvector(Kll_, Mll_, ndof_, neigenvalue)
        nmode = len(eigenvalue)
        log.debug(f"eigenvalue = {eigenvalue}")

        # Expand l-set eigenvectors back to a-set then g-set
        ndof_a = len(xa)
        phia = np.zeros((nmode, ndof_a), dtype=fdtype)
        for imode in range(nmode):
            phil = np.zeros(nl, dtype=fdtype)
            phil[ipositive] = xl_[:, imode]
            phia[imode, lset_local] = phil

        phig = np.full((nmode, ndof_work), 0.0, dtype=fdtype)
        phig[:, aset] = phia
        phig[:, sset] = xs
        out["rset_b"] = rset_b
    else:
        # No SUPORT: solve on the full a-set (original path)
        Kaa_, ipositive, unused_inegative, unused_sz_set = remove_rows(
            Kaa, aset)
        Maa_ = Maa[ipositive, :][:, ipositive]
        ndof_ = Kaa_.shape[0]

        if issparse(Maa_):
            Maa_ = Maa_.toarray()
        eigenvalue, xa_ = solve_eigenvector(Kaa_, Maa_, ndof_, neigenvalue)
        nmode = len(eigenvalue)
        log.debug(f"eigenvalue = {eigenvalue}")

        ndof_a = len(xa)
        phia = np.zeros((nmode, ndof_a), dtype=fdtype)
        for imode in range(nmode):
            phia[imode, ipositive] = xa_[:, imode]

        phig = np.full((nmode, ndof_work), 0.0, dtype=fdtype)
        phig[:, aset] = phia
        phig[:, sset] = xs

    # ----------------------MPC recovery (n -> g) for modes----------------------
    if GMN is not None:
        # phig is currently in n-set (ndof_work = n_ndof), expand to g-set
        GMN_dense = GMN.toarray()
        phig_g = np.zeros((nmode, ndof), dtype=fdtype)
        for imode in range(nmode):
            phig_g[imode, :] = (GMN_dense @ phig[imode, :]).ravel()
        phig = phig_g

    # phit = phig
    nnode_g = len(node_gridtype)
    phi, Mhh, Khh = apply_phi_normalization(Mgg, Kgg, eigenvalue, phig, nmode, nnode_g, norm_str)
    log.info(f"Mhhp_diag: {np.diag(Mhh)}")
    log.info(f"Khh_diag: {np.diag(Khh)}")

    assert np.all(np.isfinite(phia))
    assert np.all(np.isfinite(phig))
    out["modes_eigenvalue"] = eigenvalue
    out["modes_norm"] = norm_str
    out["modes_phig"] = phig
    out["modes_Mhh"] = Mhh
    out["modes_Khh"] = Khh
    return out


def _write_mass_participation_f06(f06_file: TextIO,
                                  mpf: dict[str, np.ndarray],
                                  cycle: np.ndarray,
                                  nmode: int) -> None:
    """Write modal effective mass / participation factor table to f06."""
    directions = ['T1', 'T2', 'T3', 'R1', 'R2', 'R3']
    eff_mass = mpf['effective_mass']
    eff_ratio = mpf['effective_mass_ratio']
    cum_ratio = mpf['cumulative_ratio']

    f06_file.write('\n')
    f06_file.write('                              M O D A L   E F F E C T I V E   M A S S\n')
    f06_file.write('\n')
    f06_file.write('  MODE    FREQUENCY        T1             T2             T3'
                   '             R1             R2             R3\n')

    for i in range(nmode):
        f06_file.write(f'  {i+1:4d}  {cycle[i]:12.6f}')
        for j in range(6):
            f06_file.write(f'  {eff_mass[i, j]:13.6E}')
        f06_file.write('\n')

    f06_file.write('\n')
    f06_file.write(
        '                      M O D A L   E F F E C T I V E   M A S S'
        '   F R A C T I O N\n'
    )
    f06_file.write('\n')
    f06_file.write('  MODE    FREQUENCY        T1             T2             T3'
                   '             R1             R2             R3\n')

    for i in range(nmode):
        f06_file.write(f'  {i+1:4d}  {cycle[i]:12.6f}')
        for j in range(6):
            f06_file.write(f'  {eff_ratio[i, j]:13.6E}')
        f06_file.write('\n')

    f06_file.write('\n  CUMULATIVE:\n')
    f06_file.write(f'  {"":4s}  {"":12s}')
    for j in range(6):
        f06_file.write(f'  {cum_ratio[-1, j]:13.6E}')
    f06_file.write('\n\n')

    total_mass = mpf['total_mass']
    f06_file.write('  TOTAL MASS:       ')
    for j in range(6):
        f06_file.write(f'  {total_mass[j]:13.6E}')
    f06_file.write('\n\n')


def _write_gpforce_balance(
    f06_file: TextIO,
    op2: OP2,
    node_gridtype: np.ndarray,
    Fg_applied: np.ndarray,
    fspc: np.ndarray,
    Kgg: np.ndarray,
    xg: np.ndarray,
    isubcase: int,
    ngrid: int,
    ndof_per_grid: int,
    title: str = '',
    subtitle: str = '',
    label: str = '',
    page_stamp: str = 'PAGE %s',
    page_num: int = 1,
) -> int:
    """Build grid point force balance and write to F06 via pyNastran object.

    For each grid, reports: applied load, SPC force, and the total
    (should be ~0 at equilibrium).
    """
    from pyNastran.op2.tables.ogf_gridPointForces.ogf_objects import (
        RealGridPointForcesArray,
    )

    # Internal forces: F_int = K @ x (at each DOF)
    if hasattr(Kgg, 'dot'):
        f_internal = Kgg.dot(xg)
    else:
        f_internal = Kgg @ xg

    # Clean up NaN in applied forces (SPCD DOFs)
    Fg_clean = np.where(np.isnan(Fg_applied), 0.0, Fg_applied)

    # Build arrays for RealGridPointForcesArray
    node_element_list = []
    element_names_list = []
    data_list = []

    for i in range(ngrid):
        nid = int(node_gridtype[i, 0])
        i0 = i * ndof_per_grid
        i1 = i0 + ndof_per_grid

        fa = Fg_clean[i0:i1]
        fs = fspc[i0:i1]
        fi = f_internal[i0:i1]

        max_force = max(
            np.max(np.abs(fa)),
            np.max(np.abs(fs)),
            np.max(np.abs(fi)))
        if max_force < 1e-20:
            continue

        if np.max(np.abs(fa)) > 1e-20:
            node_element_list.append([nid, 0])
            element_names_list.append('APP-LOAD')
            data_list.append(fa[:6])

        if np.max(np.abs(fs)) > 1e-20:
            node_element_list.append([nid, 0])
            element_names_list.append('F-OF-SPC')
            data_list.append(fs[:6])

        if np.max(np.abs(fi)) > 1e-20:
            node_element_list.append([nid, 0])
            element_names_list.append('F-OF-SPC')
            data_list.append(fi[:6])

        total = fa + fs - fi
        node_element_list.append([nid, 0])
        element_names_list.append('*TOTALS*')
        data_list.append(total[:6])

    if not node_element_list:
        return page_num

    node_element = np.array(node_element_list, dtype='int32')
    element_names = np.array(element_names_list, dtype='U8')
    data = np.array(data_list, dtype='float32').reshape(1, -1, 6)

    gpf_obj = RealGridPointForcesArray.add_static_case(
        'OGPFB1', node_element, element_names, data,
        isubcase, title=title, subtitle=subtitle, label=label)
    op2.grid_point_forces[isubcase] = gpf_obj

    page_num = gpf_obj.write_f06(
        f06_file, header=None, page_stamp=page_stamp,
        page_num=page_num, is_mag_phase=False, is_sort1=True)
    return page_num


def solve_eigenvector(
    Kaa: csc_matrix,
    Maa: np.ndarray,
    ndof: int,
    neigenvalues: int,
    use_lobpcg: bool = False,
    X0: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray]:
    """Solve the generalized eigenproblem K*x = lambda*M*x.

    Parameters
    ----------
    use_lobpcg : bool
        If True, use LOBPCG instead of ARPACK (eigsh). Requires pyamg
        for preconditioning (optional but recommended).
    X0 : ndarray, optional
        Initial eigenvector guess for LOBPCG warm-start, shape (ndof, k).
        Pass previous converged modes for fast re-solves.
    """
    Kaa_dense = Kaa.toarray()
    assert isinstance(Maa, np.ndarray), type(Maa)

    if not np.any(Maa):
        from pyNastran.f06.errors import FatalError
        raise FatalError('mass matrix is zero; check that elements have material density')

    # reduce the size of the matrix going into the solver
    # by removing empty rows/columns
    Kaai = Kaa_dense
    is_kaa = np.abs(Kaai) > 0.0
    is_maa = np.abs(Maa) > 0.0
    is_kaa_maa = is_kaa & is_maa
    rows = np.any(is_kaa_maa, axis=0)
    cols = np.any(is_kaa_maa, axis=1)
    assert len(rows) == len(cols)
    is_modes = rows & cols
    no_modes = ~is_modes
    Kaa2 = Kaa[is_modes, :][:, is_modes]
    Maa2 = Maa[is_modes, :][:, is_modes]

    ndof2 = Maa.shape[0]
    backend = get_solver()
    if ndof2 <= neigenvalues:
        Kaa2_dense = Kaa2.todense()
        eigenvalues, xa = sp.linalg.eigh(Kaa2_dense, Maa2)
    elif use_lobpcg:
        X0_reduced = None
        if X0 is not None:
            X0_reduced = X0[is_modes, :]
        eigenvalues, xa_ = backend.lobpcg(Kaa2, k=neigenvalues, M=Maa2, X0=X0_reduced)
        xa = np.full((ndof, neigenvalues), np.nan, dtype=xa_.dtype)
        xa[no_modes, :] = 0
        xa[is_modes, :] = xa_
    else:
        from scipy.sparse.linalg import ArpackNoConvergence
        try:
            eigenvalues, xa_ = backend.eigsh(
                Kaa2, k=neigenvalues, M=Maa2, which="SM", return_eigenvectors=True
            )
        except ArpackNoConvergence as e:
            from pyNastran.f06.errors import FatalError
            raise FatalError(
                f'eigensolver did not converge: {e.args[0]}'
            ) from e
        xa = np.full((ndof, neigenvalues), np.nan, dtype=xa_.dtype)
        xa[no_modes, :] = 0
        xa[is_modes, :] = xa_
    return eigenvalues, xa


def _build_xg(model: BDF, dof_map: DOF_MAP, ndof: int, subcase: Subcase) -> NDArrayNfloat:
    """
    Builds the {xg} vector, which has all SPCs in the analysis (cd) frame
    (called global g by NASTRAN)

    {s} = {sb} + {sg}
    {sb} = SPC set on SPC/SPC1/SPCADD cards (boundary)
    {sg} = SPCs from PS field on GRID card (grid)

    """
    log = model.log
    # model = self.model
    xspc = np.full(ndof, np.nan, dtype="float64")
    if "SPC" not in subcase:
        log.warning("no SPCs defined in case control")
        spc_set = np.array([], dtype="int32")
        sset = np.zeros(ndof, dtype="bool")
        return spc_set, sset, xspc
    spc_id, unused_options = subcase["SPC"]
    spcs = []
    # for spc in model.spcs:
    # if spc.n == 0:
    # continue
    # spci = spc.slice_card_by_id(spc_id)
    # spcs.append(spci)
    spc_cards = [spc for spc in model.spc_cards if spc.n > 0]
    spcs = [spc.slice_card_by_id(spc_id, sort_ids=True) for spc in spc_cards]
    model.spc1
    # spcs = model.get_reduced_spcs(spc_id, consider_spcadd=True, stop_on_failure=True)

    spc_set = []
    sset = np.zeros(ndof, dtype="bool")
    for spc in spcs:
        if spc.type == "SPC1":
            # print(spc.get_stats())
            dofs_missed = []
            for comp, (inode0, inode1) in zip(spc.components, spc.inode):
                for dof in str(comp):
                    dofi = int(dof)
                    spc_nodes = spc.node_id[inode0:inode1]
                    for nid in spc_nodes:
                        try:
                            idof = dof_map[(nid, dofi)]
                        except Exception:
                            dofs_missed.append((nid, dofi))
                            # print('dof_map =', dof_map)
                            # print((nid, dofi))
                            continue
                        sset[idof] = True
                        spc_set.append(idof)
                        xspc[idof] = 0.0
            if dofs_missed:
                dof_str = ", ".join(str(dofi) for dofi in dofs_missed)
                log.warning(f"Missing (nid,dof) pairs:")
                log.warning(f"  {dof_str}\n{spc.rstrip()}")
                # =({nid},{dofi}) doesn't exist...skipping
            del dofs_missed
        elif spc.type == "SPC":
            for nid, components, enforcedi in zip(spc.node_id, spc.components, spc.enforced):
                for component in str(components):
                    dofi = int(component)
                    idof = dof_map[(nid, dofi)]
                    sset[idof] = True
                    spc_set.append(idof)
                    xspc[idof] = enforcedi
        else:
            raise NotImplementedError(spc)
    spc_set = np.array(spc_set, dtype="int32")
    # print('spc_set =', spc_set, xspc)
    return spc_set, sset, xspc
