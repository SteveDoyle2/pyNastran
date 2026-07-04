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
import scipy
from scipy.linalg import eigh, eig
from scipy.sparse import csc_matrix, lil_matrix, dok_matrix, issparse
from scipy.sparse.linalg import ArpackNoConvergence

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
#from pyNastran.utils.numpy_utils import integer_types  # , float_types
from pyNastran.dev.bdf_vectorized3.bdf import BDF, Subcase
from pyNastran.dev.bdf_vectorized3.cards.constraints import MPC

from pyNastran.f06.f06_writer import make_end
from pyNastran.f06.f06_tables.oload_resultant import Resultant
from pyNastran.f06.errors import FatalError

from pyNastran.op2.op2 import OP2
from pyNastran.op2.op2_interface.op2_classes import (
    RealDisplacementArray,
    RealSPCForcesArray,
    RealLoadVectorArray,
    ComplexDisplacementArray,
    ComplexVelocityArray,
    ComplexAccelerationArray,
    RealEigenvectorArray,
    RealMPCForcesArray,
    RealGridPointForcesArray,
)
from pyNastran.op2.result_objects.grid_point_weight import make_grid_point_weight
# from pyNastran.bdf.mesh_utils.loads import get_ndof

from .loads.reduce_pg import reduce_Pg_to_Pa
from .recover.table import (
    _save_displacment, _save_spc_forces,
    _save_mpc_forces, _save_applied_load)

from .recover.freq_force import recover_force_freq
from .recover.modal_force import recover_force_103
from .recover.static_force import recover_force_101

from .recover.static_stress import recover_stress_101
from .recover.static_strain import recover_strain_101
from .recover.strain_energy import recover_strain_energy_101
from .recover.utils import (
    get_f06_op2_pch_set, get_mag_phase_from_options, get_plot_request)

from .partition import (
    partition_matrix, # partition_vector,
    partition_vector2, partition_vector3)
from .utils_modes import (
    slice_modal_set, get_real_eigenvalue_method,
    apply_phi_normalization,
    save_eigenvalues,
    compute_mass_participation,
)
from .modal_frequency import get_freq_damping
from .utils_freq import get_frequencies, slice_freq_set

#-----------------------------------------------
from .matrices.build_stiffness import (
    build_Kgg, Kbb_to_Kgg)
from .matrices.build_stiffness_geometric import build_KDgg
from .matrices.build_mass import build_Mbb
from .loads.build_fb import build_Fb_from_loadid
#--------------------------------------------------------
from .utils import recast_data, get_param, DOF_MAP

from .craig_bampton import (
    run_craig_bampton, write_cb_to_op4, write_cb_to_h5)

from pyNastran.dev.bdf_vectorized3.mesh_utils.inertia_relief import (
    build_rigid_body_modes, compute_inertia_relief)
from pyNastran.dev.bdf_vectorized3.solver.matrices.gmn_matrix import assemble_gmn
from .dynamic_file_writer import DynamicFileWriter

Array = np.ndarray | scipy.sparse.csc_matrix


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
        self.Mgg_override = None
        self.KDgg_override = None

        base_name = os.path.splitext(model.bdf_filename)[0]
        self._bdf_filename = base_name + ".solver.bdf"
        self.f06_filename = base_name + ".solver.f06"
        self.op2_filename = base_name + ".solver.op2"
        self.op4_filename = base_name + ".solver.op4"
        self.h5_filename = base_name + ".solver.h5"
        # print(self.f06_filename)
        # print(self.op2_filename)

        self.solver_dict = {
            'op4': DynamicFileWriter(self.op4_filename),
            'h5': DynamicFileWriter(self.h5_filename),
        }
        #print(self.solver_dict)


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
                title = subcase.get_parameter("TITLE")[0]
                assert isinstance(title, str), title
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
            else:  # pragma: no cover
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
        for name, file_obj in self.solver_dict.items():
            file_obj.close()
        self.solver_dict = {}

    def _update_card_count(self) -> None:
        for card_type, values in self.model._type_to_id_map.items():
            self.model.card_count[card_type] = len(values)

    def build_Fb(self, xg: NDArrayNfloat, sset_b,
                 dof_map: DOF_MAP, ndof: int,
                 subcase: Subcase,
                 Mbb: Array | None = None,
                 xyz_cid0: np.ndarray | None = None) -> NDArrayNfloat:
        model = self.model
        log = model.log
        log.info("starting build_Fb")

        has_load = "LOAD" in subcase
        has_temp = "TEMPERATURE(LOAD)" in subcase or "TEMPERATURE(BOTH)" in subcase

        load_id = 0
        temp_load_id = 0
        if has_load:
            load_id, unused_options = subcase["LOAD"]

        if has_temp:
            if "TEMPERATURE(LOAD)" in subcase:
                temp_load_id, _ = subcase["TEMPERATURE(LOAD)"]
            else:
                temp_load_id, _ = subcase["TEMPERATURE(BOTH)"]

        if has_load or has_temp:
            log.warning('creating Fb')

        Fb, Mbb, xyz_cid0 = build_Fb_from_loadid(
            model, dof_map, ndof,
            xg, sset_b,
            load_id=load_id, temp_load_id=temp_load_id,
            Mbb=Mbb, xyz_cid0=xyz_cid0)
        log.info("end of build_Fb")
        return Fb, Mbb, xyz_cid0

    def build_GMN(self, subcase: Subcase,
                  dof_map: DOF_MAP, ndof: int,
                  xyz_cid0: np.ndarray,
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
        mpc_id = 0
        has_mpc = "MPC" in subcase
        if has_mpc:
            mpc_id, unused_options = subcase["MPC"]
        GMN_csc, mset = _build_Gmn(
            self.model, mpc_id,
            dof_map, ndof, xyz_cid0, fdtype=fdtype,
            solver_dict=self.solver_dict)
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
        fdtype: str = "float64",) -> Any:
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
        out = {}
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
        idtype = 'int32'
        # -----------------------------------------------------------------------
        model = self.model
        log = model.log

        model.setup(run_geom_check=True)
        xyz_cid0 = model.grid.xyz_cid0()
        # -----------------------------------------------------------------------
        load_id, spc_id, mpc_id, suport_id, method_id = get_case_control(
            model.sol, subcase)
        # -----------------------------------------------------------------------

        is_mpc_request = get_plot_request(subcase, 'MPCFORCES')[-1]
        is_spc_request = get_plot_request(subcase, 'SPCFORCES')[-1]
        is_gpforce_request = get_plot_request(subcase, 'GPFORCE')[-1]
        is_oload_request = get_plot_request(subcase, 'OLOAD')[-1]

        dof_map, ps = _get_dof_map(model)

        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        ngrid, ndof_per_grid, ndof = get_ndof(model, subcase)

        gset_b = ps_to_sg_set(ndof, ps)
        rset_b = get_rset_bool(model, dof_map, ndof, suport_id=suport_id)
        has_suport = np.any(rset_b)
        #--------------------------------------------------
        log.warning('creating Kgg')
        Kgg = get_Kgg(
            model, dof_map, ndof, ngrid, ndof_per_grid,
            idtype=idtype, fdtype=fdtype,
            Kgg_override=self.Kgg_override,
            solver_dict=self.solver_dict)

        Mbb = None
        is_mass_load = len(model.grav) > 0 or len(model.rforce) > 0
        is_param_grdpnt = 'GRDPNT' in model.params

        if is_mass_load or is_param_grdpnt or has_suport:
            log.warning('creating Mbb')
            Mbb = build_Mbb(
                model, subcase, dof_map, ndof, fdtype=fdtype,
                solver_dict=self.solver_dict)
            page_num = write_grid_point_weight(
                model, Mbb, dof_map, ndof,
                self.op2, f06_file,
                title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp, page_num=page_num)

        sset, sset_b, xg = _build_xg(model, dof_map, ndof, subcase)

        Fb, Mbb, xyz_cid0 = self.build_Fb(
            xg, sset_b, dof_map, ndof, subcase, Mbb=Mbb,
            xyz_cid0=xyz_cid0)

        # ---------------------MPC reduction (g -> n)---------------------
        GMN, mset = self.build_GMN(
            subcase, dof_map, ndof, xyz_cid0, fdtype=fdtype)
        self.GMN = GMN
        self.mset = mset
        #out['GMN'] = GMN
        out['mset_bool'] = mset

        log.warning('creating Mgg')
        Mgg = None if Mbb is None else Kbb_to_Kgg(
            model, Mbb, ngrid, ndof_per_grid,
            inplace=False)
        write_mat(model, Mgg, 'PRTMGG', self.solver_dict)

        Knn, Mnn, autospc_n_dofs = static_g_to_n(
            model,
            Kgg, Mgg, GMN,
            log, self.solver_dict)
        self._autospc_n_dofs = autospc_n_dofs
        gset = np.arange(ndof, dtype=idtype)

        write_mat(model, Fb, 'PRTPB', self.solver_dict)
        Fg = xb_to_xg(model, Fb, ngrid, ndof_per_grid)
        write_mat(model, Fg, 'PRTPG', self.solver_dict)

        # Save g-set force for oload output before MPC transform
        Fg_gset = np.where(np.isnan(Fg), 0.0, Fg.copy())

        asetmap = get_aset(model)

        Fn, nset_indices = static_Pg_to_Pn(
            model,
            Fg, GMN,
            gset, mset,
            log, self.solver_dict)

        if GMN is not None:
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

            if asetmap:
                aset_g = apply_dof_map_to_set(asetmap, dof_map, idtype=idtype)
                aset = np.array([
                    g_to_n[i] for i in aset_g if g_to_n[i] >= 0],
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

            # -------------------SPC reduction (g -> a)-------------------
            if asetmap:
                aset = apply_dof_map_to_set(
                    asetmap, dof_map, idtype=idtype)
            else:
                aset = np.setdiff1d(gset, sset)

        #if not is_gpforce_request:
        #    del Kgg_orig

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
        page_num = write_oload_resultant(
            Fg_oload, dof_map, isubcase, ngrid, ndof_per_grid, f06_file, page_stamp, page_num, log)

        # aset - analysis set
        # sset - SPC set
        # print('aset = ', aset)
        # print('sset = ', sset)

        # u1 = Kaa^-1 * (F1k - Kas*u2k)
        abs_xg = np.abs(xg)
        finite_xg = np.any(np.isfinite(abs_xg))
        if finite_xg and np.nanmax(abs_xg) > 0.0:
            self.log.warning(f"SPCD found")
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

        log.warning('partitioning F: a = g - s - 0')
        Fa, Fs = partition_vector2(Fg, [["a", aset], ["s", sset]])
        del Fs

        log.warning('partitioning U: a = g - s - 0')
        xa, xs, x0 = partition_vector3(xg, [["a", aset], ["s", sset], ["0", set0]])
        write_mat(model, xs, 'PRTUS', self.solver_dict)
        write_mat(model, x0, 'PRTU0', self.solver_dict)
        # write_mat(model, xa, 'PRTUA', self.solver_dict)
        # self.log.info(f'xg = {xg}')
        del xg
        # self.log.info(f'xa = {xa}')
        # self.log.info(f'xs = {xs}')
        # self.log.info(f'x0 = {x0}')
        del x0

        self.Kgg = Kgg_orig
        log.warning('partitioning K: a = g - s - 0')
        K = partition_matrix(Kgg, [("a", aset), ("s", sset), ("0", set0)])
        Kaa = K["aa"]
        Kss = K["ss"]
        # Kas = K['as']
        Ksa = K["sa"]
        write_mat(model, Kaa, 'PRTKAA', self.solver_dict)
        write_mat(model, Kss, 'PRTKSS', self.solver_dict)
        write_mat(model, Ksa, 'PRTKSA', self.solver_dict)
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

        write_mat(model, Fa_solve, 'PRTFA2', self.solver_dict)

        # --- SUPORT / inertia relief for statics ---
        inrel = get_param(model, 'INREL', -1)

        if has_suport and inrel == -2 and is_aset:
            # Inertia relief: partition a = l + r, apply inertia relief
            log.warning('SUPORT/INREL: l = a - r')
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
            grid_ids = model.grid.node_id
            D_full = build_rigid_body_modes(grid_ids, xyz_cid0)

            # Extract D for a-set DOFs only
            D_a = D_full[a_indices, :]
            write_mat(model, D_a, 'PRTDRB', self.solver_dict)

            # Build a-set mass matrix (dense)
            M_a = partition_matrix(
                Mgg, [("a", aset), ("s", sset), ("0", set0)])
            Maa = M_a["aa"]
            write_mat(model, Maa, 'PRTMAA', self.solver_dict)
            Maa_dense = todense(Maa)
            del Maa
                
            # Apply inertia relief
            F_net, a_rigid, F_inertia = compute_inertia_relief(
                Maa_dense, D_a, Fa_solve)
            log.info(f"  rigid body acceleration: {a_rigid}")

            # Partition to l-set and solve
            Kaa_dense = todense(Kaa)
            Kll = Kaa_dense[np.ix_(lset_local, lset_local)]
            Fl = F_net[lset_local]
            write_mat(model, Kll, 'PRTKLL', self.solver_dict)
            write_mat(model, Fl, 'PRTPL', self.solver_dict)

            # Solve Kll @ ul = Fl
            Kll_sp = csc_matrix(Kll)
            xl_, ipositive_l, inegative_l = solve(
                Kll_sp, Fl, np.ones(nl, dtype='bool'), log,
                idtype=idtype)
            write_mat(model, xl_, 'PRTUL', self.solver_dict)

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
            log.warning('a set reduction')
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
            FatalError("A-set is empty; all DOFs are constrained")
            #self.xa_ = []
            #self.Fa_ = []

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
            log.debug(f"  fspc_0 = K0a @ xa + K0s @ xs")
            fspc_0 = K0a @ xa + K0s @ xs
            log.debug(f"  fspc_0 = {fspc_0}")
            Fg[set0] = fspc_0
            fspc[set0] = fspc_0

        # ---------------------MPC recovery (n -> g)---------------------
        #is_load_recovery_request = (
        #    is_spc_request or is_mpc_request or
        #    is_gpforce_request # or is_oload_request
        #    or 1
        #)

        if GMN is not None:
            log.warning('recover SPC/MPC forces')
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
        #write_mat(model, xs, 'PRTUB', self.solver_dict)
        #write_mat(model, xs, 'PRTPB', self.solver_dict)

        # log.debug(f'Fs = {Fs}')
        # log.debug(f'Fb = {Fb}')
        # log.debug(f'xb = {xb}')

        # SPCFORCE resultant
        #if 'SPCFORCES' in subcase:
        ndofi = ngrid * ndof_per_grid
        page_num = _write_spcforce_resultant(
            f06_file, fspc, ngrid, ndof_per_grid, isubcase,
            page_stamp, page_num, log)

        op2 = self.op2
        page_stamp % page_num
        page_num = recover_statics(
            model, op2, f06_file,
            subcase,
            xg,
            Fg_oload, fspc, self.fmpc,
            node_gridtype,
            ngrid, ndof_per_grid,
            title=title, subtitle=subtitle, label=label,
            page_stamp=page_stamp, page_num=page_num,
            fdtype=fdtype)

        if 'GPFORCE' in subcase:
            page_num = _write_gpforce_balance(
                f06_file, self.op2, node_gridtype, Fg_oload, fspc,
                Kgg_orig, xg, isubcase, ngrid, ndof_per_grid,
                title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp, page_num=page_num)

        page_stamp += "\n"
        page_num = recover_force_101(
            f06_file,
            op2,
            self.model,
            dof_map,
            subcase,
            xb,
            title=title,
            subtitle=subtitle,
            label=label,
            page_stamp=page_stamp)

        page_num = recover_strain_101(
            f06_file, op2, self.model, dof_map, subcase, xb,
            title=title, subtitle=subtitle, label=label,
            page_stamp=page_stamp)
        page_num = recover_stress_101(
            f06_file, op2, self.model, dof_map, subcase, xb,
            title=title, subtitle=subtitle, label=label,
            page_stamp=page_stamp)
        page_num = recover_strain_energy_101(
            f06_file,
            op2,
            self.model,
            dof_map,
            subcase,
            xb,
            title=title,
            subtitle=subtitle,
            label=label,
            page_stamp=page_stamp,)
        self.log.info("finished")
        out = {}
        return out, page_num, end_options

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
        matrices_to_save: set[str] | None=None):
        """
        [M]{xdd} + [C]{xd} + [K]{x} = {F}
        [M]{xdd} + [K]{x} = {0}
        -[M]{xdd}λ^2 + [K]{x} = {0}
        {X}(λ^2 - [M]^-1[K]) = {0}
        λ^2 - [M]^-1[K] = {0}
        λ^2 = [M]^-1[K]
        [A][X] = [X]λ^2
        """
        if matrices_to_save is None:
            matrices_to_save = set([])

        model = self.model
        nmodes, norm_str = get_real_eigenvalue_method(model, subcase)
        log = model.log
        log.debug("run_sol_103 (modes)")
        assert len(model.methods), "SOL 103 (modes) requires a METHOD and a EIGR/EIGRL card"
        end_options = [
            "SEMR",  # MASS MATRIX REDUCTION STEP (INCLUDES EIGENVALUE SOLUTION FOR MODES)
            "SEKR",  # STIFFNESS MATRIX REDUCTION STEP
            "MODES",  # run modes
        ]
        op2 = self.op2
        # write_f06 = True

        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        xyz_cid0 = model.grid.xyz_cid0()
        # -----------------------------------------------------------------------
        load_id, spc_id, mpc_id, suport_id, method_id = get_case_control(
            model.sol, subcase)
        # -----------------------------------------------------------------------

        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof = get_ndof(self.model, subcase)

        gset_b = ps_to_sg_set(ndof, ps)
        rset_b = get_rset_bool(model, dof_map, ndof, suport_id=suport_id)
        has_suport = np.any(rset_b)
        #--------------------------------------------------

        # Build GMN for MPC reduction
        GMN, mset = self.build_GMN(
            subcase, dof_map, ndof, xyz_cid0, fdtype=fdtype)
        self.GMN = GMN
        self.mset = mset

        matrices_to_save.add('Mgg')
        page_num, out = _run_modes(
            model, subcase,
            op2, f06_file,
            ngrid,
            ndof_per_grid,
            ndof,
            node_gridtype,
            dof_map,
            Kgg_override=self.Kgg_override,
            Mgg_override=self.Mgg_override,
            GMN=GMN,
            mset=mset,
            idtype=idtype,
            fdtype=fdtype,
            page_stamp=page_stamp, page_num=page_num,
            solver_dict=self.solver_dict,
            matrices_to_save=matrices_to_save,
        )
        if 'GMN' in matrices_to_save:
            out['GMN'] = GMN
            out['mset_bool'] = mset

        phig = out["modes_phig"]
        eigenvalue = out["modes_eigenvalue"]
        cycle = np.sqrt(np.abs(eigenvalue)) / (2.0 * np.pi)
        mode_cycle = eigenvalue
        isubcase = subcase.id

        page_num = save_eigenvalues(
            self.op2, f06_file, out, subcase, title,
            page_stamp=page_stamp, page_num=page_num)
        eigenvalue = out["modes_eigenvalue"]
        nmode = len(eigenvalue)

        # Mass participation factors
        Mgg = out["Mgg"]
        nnode_g = node_gridtype.shape[0]
        mpf = compute_mass_participation(phig, Mgg, nnode_g, nmode)
        out["mass_participation"] = mpf

        _write_mass_participation_f06(f06_file, mpf, cycle, nmode)
        log.info(f"mass participation cumulative (Tx): {mpf['cumulative_ratio'][-1, 0]:.4f}")
        page_num = _save_eigenvectors(
            model, subcase,
            phig, eigenvalues, mode_cycle, node_gridtype,
            op2, f06_file, page_num, sol_type='modes',
            page_stamp=page_stamp,
            title=title, subtitle=subtitle, label=label,)

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

        # add transform for xg to xb
        # (need for rods, quads, etc., but not springs)
        phib = xg_to_xb(model, phig, ngrid, ndof_per_grid)
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

        # Per-mode force/stress/strain/ESE recovery
        # using static routines
        #
        # TODO: how does this write correctly to the op2?
        # TODO: this does a lot of extra work 
        nmode = phig.shape[0]
        for imode in range(nmode):
            xb_mode = phib[imode, :]
            recover_force_101(
                f06_file, op2, self.model, dof_map, subcase,
                xb_mode, title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp,)
            recover_strain_101(
                f06_file, op2, self.model, dof_map, subcase,
                xb_mode, title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp,)
            recover_stress_101(
                f06_file, op2, self.model, dof_map, subcase,
                xb_mode, title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp,)
            recover_strain_energy_101(
                f06_file, op2, self.model, dof_map, subcase,
                xb_mode, title=title, subtitle=subtitle, label=label,
                page_stamp=page_stamp,)
        str(f06_file)
        str(page_stamp)
        return out, page_num, end_options

    def build_complex_Fb(self,
                         model: BDF,
                         subcase: Subcase,
                         ndof_g: int):
        # TODO: parse inputs for loads
        Fg = np.ones((ndof_g, 1), dtype="complex64")
        dload_id, _ = subcase['DLOAD']

        return Fg

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
        fdtype: str = "float64"):
        model = self.model
        neigenvalue, norm_str = get_real_eigenvalue_method(model, subcase)
        if norm_str == "MAX":
            raise RuntimeError("norm_str=MAX and should be MASS (it makes the math harder)")

        log = model.log
        isubcase = subcase.id
        op2 = self.op2
        # ---------------------------------------------------------------------
        # case control
        load_id, spc_id, mpc_id, suport_id, method_id = get_case_control(
            model.sol, subcase)
        dload_id, _ = subcase['DLOAD']
        freq_id, _ = subcase['FREQUENCY']
        # ---------------------------------------------------------------------
        ngrid = len(model.grid)
        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        ngrid, ndof_per_grid, ndof_g = get_ndof(self.model, subcase)

        # handles MAX/MASS normalization
        log.info('running modes')
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
            matrices_to_save={'Koa', 'Koo', 'GMN'},
        )
        Koa = out.get('Koa')
        Koo = out.get('Koo')
        
        GMN = None
        if mpc_id or 'GMN' in out:
            GMN = out.get('GMN')
            mset_b = out['mset_bool']
        sset_b = out['sset_bool']
        oset_b = None

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
        assert eigenvalue.ndim == 1, eigenvalue
        omegan = np.sqrt(np.abs(eigenvalue))

        # frequency
        log.info('getting frequencies')
        freq = get_frequencies(model, subcase, omegan)
        nfreq = len(freq)
        omega = 2 * np.pi * freq
        assert nfreq > 0, nfreq

        Fb = self.build_complex_Fb(model, subcase, ndof_g)
        Fg = xb_to_xg(model, Fb, ngrid, ndof_per_grid)
        Fa = reduce_Pg_to_Pa(
            Fg,
            ndof_g,
            sset_b,
            mset_b=mset_b,
            oset_b=oset_b,
            Gm=GMN, Koo=Koo, Koa=Koa)

        Cstr_gg, zomegan2 = get_freq_damping(model, omegan, ndof_g)

        # nfreq = len(omegas)
        # shape = (ndof, ndof, nfreq)

        # phiT = phi.T
        #  modal space (h); sometimes called x
        # print(phi.shape, Mgg.shape)
        # phi = phit.T
        Mhh = out["modes_Mhh"]
        Khh = out["modes_Khh"]
        assert Khh.ndim == 2, Khh.shape
        assert Mhh.ndim == 2, Mhh.shape
        # Mhh = phit @ Mgg @ phi
        # Chh = phi @ Cgg @ phi.T
        # Khh = phi @ Kgg @ phi.T
        # Fh = phi @ Fg

        # U(t)   = A*e^(i*omega*t) U
        # Ud(t)  = i*omega * A*e^(i*omega*t) U
        # Udd(t) = -omega^2 * A*e^(i*omega*t) U

        #log.info('[Mgg] {xddg} + Cgg xdg + Kgg xg = Fg')
        #log.info('{qh}   = [PHI]{xg}')
        #log.info('{qdh}  = [PHI]{xdg}')
        #log.info('{qddh} = [PHI]{xddg}')
        #log.info('[PHIT][Mgg][PHI]{qddh} + [PHIT][Cgg][PHI]{qdg} + '
        #         '[PHIT][Kgg][PHI]{qh} = [PHIT]{Fg}')
        #log.info('[I]{qddh} + [C]{qdh} + [omega]{qh} = {Fh}')
        
        Fh = phit @ Fg
        # print('Fh:\n', Fh)
        omega2 = omega * omega  # frequncy
        omegan2 = np.diag(Khh)  # natural frequency squared
        omegan = np.sqrt(omegan2)

        khh_diag = np.diag(Khh)  # TODO: this part is weird
        p = 1.0
        # print(f'ndof_g={ndof_g} phi.shape={phi.shape}')
        xq = np.zeros((nfreq, nmode), dtype="complex64")
        for ifreq, omega2i in enumerate(omega2):
            log.info(f"ifreq={ifreq}")
            # denom = np.sqrt((1-omega2i/omegan2)**2 + zomegan2**2)
            # p / (k * denom) * np.sin(omegan)

            # $$ u(t) =
            tf = np.sqrt((1 - omega2i / omegan2) ** 2 + (zomegan2 / omegan) ** 2)
            mag = p / khh_diag * tf
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
                log.info(f"skipping {key!r} because no F06/OP2 output was requested")
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
            if write_op2:
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
                    is_sort1=True,)
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
        write_op2: bool = True,):
        """SOL 105 linear buckling: (K + lambda*KD)*x = 0.

        Subcase must have LOAD, SPC, and METHOD.
        STATSUB(PRELOAD) is not yet supported — uses the same subcase for preload.
        """
        model = self.model
        log = model.log
        op2 = self.op2

        log.debug("run_sol_105 (buckling)")
        neigenvalue, norm_str = get_real_eigenvalue_method(
            model, subcase)

        end_options = ["SEKR", "MODES"]

        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof = get_ndof(model, subcase)
        xyz_cid0 = model.grid.xyz_cid0()
        # -----------------------------------------------------------------------

        GMN = None
        Kgg = get_Kgg(
            model, dof_map, ndof, ngrid, ndof_per_grid,
            idtype=idtype, fdtype=fdtype,
            Kgg_override=self.Kgg_override,
            solver_dict=self.solver_dict)

        # Static preload solve
        sset, sset_b, xg = _build_xg(model, dof_map, ndof, subcase)

        Fb, Mbb, xyz_cid0 = self.build_Fb(
            xg, sset_b, dof_map, ndof, subcase,
            Mbb=None, xyz_cid0=xyz_cid0)
        Mbb = None
        Fg = xb_to_xg(model, Fb, ngrid, ndof_per_grid)

        free_dofs = np.where(~sset_b)[0]
        assert len(sset_b) == len(Fg)
        Kff = Kgg.tocsc()[np.ix_(free_dofs, free_dofs)].toarray()
        Ff = Fb[free_dofs]

        try:
            u_free = np.linalg.solve(Kff, Ff)
        except:
            print(f'Kff:\n{Kff}')
            raise
        u_global = np.zeros(ndof, dtype=fdtype)
        u_global[free_dofs] = u_free

        log.debug("  static preload solved")
        out = _run_buckling(
            model, subcase,
            ndof, dof_map, u_global, Kff, free_dofs)

        # "buckling_eigenvalues": eigenvalue,
        # "buckling_modes": xa_,
        # "free_dofs": free_dofs,
        phi_f = out['buckling_modes']
        eigenvalues = out['buckling_eigenvalues']
        
        recover_results = get_recover_results(subcase, ['DISPLACEMENT'])

        if recover_results:
            nmode = phi_f.shape[0]
            phi_n = np.zeros((nmode, ndof), dtype=fdtype)
            phi_n[:, free_dofs] = phi_f
            phi_g = phi_n_to_g(phi_n, ndof, GMN=GMN, fdtype=fdtype)

            page_num = _save_eigenvectors(
                model, subcase,
                phi_g, eigenvalues, eigenvalues, node_gridtype,
                op2, f06_file, page_num, sol_type='buckling',
                page_stamp=page_stamp,
                title=title, subtitle=subtitle, label=label,)

        self.buckling_eigenvalues = out['buckling_eigenvalues']
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
        fdtype: str = "float64",):
        """SOL 31 Craig-Bampton (GEN CB MODEL).

        Boundary DOFs defined by SUPORT/SUPORT1 cards.
        Number of modes from METHOD (EIGRL/EIGR).
        """
        model = self.model
        neigenvalues, _ = get_real_eigenvalue_method(model, subcase)

        base_name = os.path.splitext(model.bdf_filename)[0]
        op4_filename = base_name + ".op4"
        h5_filename = base_name + ".h5"

        log = model.log
        log.debug("run_sol_31 (Craig-Bampton)")

        end_options = ["SEMG", "SEMR", "SEKR"]

        result = {}
        model.setup(run_geom_check=True)
        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof = get_ndof(model, subcase)
        #-----------------------------------------------------------

        # Build global matrices
        Kgg = get_Kgg(
            model, dof_map, ndof, ngrid, ndof_per_grid,
            idtype=idtype, fdtype=fdtype,
            Kgg_override=self.Kgg_override,
            solver_dict=self.solver_dict)

        Mbb = build_Mbb(
            model, subcase, dof_map, ndof, fdtype=fdtype,
            solver_dict=self.solver_dict)
        page_num = write_grid_point_weight(
            model, Mbb, dof_map, ndof,
            self.op2, f06_file,
            title=title, subtitle=subtitle, label=label,
            page_stamp=page_stamp, page_num=page_num)
        Mgg = Kbb_to_Kgg(model, Mbb, ngrid, ndof_per_grid,
                         inplace=False)

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
        r_set_dofs = get_cb_rset_dofs(model, subcase, dof_map, g_to_f)
        
        # f - r = l?

        # Build DOF map for the free set
        f_dof_map: DOF_MAP = {}
        for (nid, dof), g_idx in dof_map.items():
            f_idx = g_to_f[g_idx]
            if f_idx >= 0:
                f_dof_map[(nid, dof)] = f_idx

        # Number of eigenvalues
        log.info(f"  R-set DOFs: {len(r_set_dofs)}, modes requested: {neigenvalues}")

        #f = n - s        unconstrained (free structural DOFs)
        #a = f - o        ananlysis??? set
        #L = a - r
        
        # Run Craig-Bampton
        cb_result = run_craig_bampton(
            Kff, Mff, f_dof_map, r_set_dofs,
            neigenvalues=neigenvalues, log=log)
        #result['Mgg'] = Mgg
        #result['Mbb'] = Mbb
        #result['Kgg'] = Kgg
        #result['Kbb'] = Kbb

        # Write summary to F06
        eigenvalues = cb_result["eigenvalues"]
        nvec = len(eigenvalues)
        f06_file.write(f"\n{'='*72}\n")
        f06_file.write("  CRAIG-BAMPTON MODEL GENERATION (SOL 31)\n")
        f06_file.write(f"{'='*72}\n\n")
        f06_file.write(f"  Number of boundary (R-set) DOFs: {len(r_set_dofs)}\n")
        f06_file.write(f"  Number of fixed-interface modes: {nvec}\n")
        f06_file.write(f"  Total CB DOFs: {cb_result['num_cb_dofs']}\n\n")

        f06_file.write("  FIXED-INTERFACE NATURAL FREQUENCIES\n")
        f06_file.write(f"  {'MODE':>4s}  {'EIGENVALUE':>14s}  {'FREQUENCY (Hz)':>14s}\n")
        for i, ev in enumerate(eigenvalues):
            freq_hz = np.sqrt(abs(ev)) / (2.0 * np.pi)
            f06_file.write(f"  {i+1:4d}  {ev:14.6e}  {freq_hz:14.6f}\n")
        f06_file.write("\n")

        # Store results on solver object
        self.cb_result = cb_result

        # Write CB matrices to OP4 and HDF5
        write_cb_to_op4(
            op4_filename, cb_result,
            is_binary=True, precision="double")
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
        fdtype: str = "float64"):
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
        neigenvalues, _ = get_real_eigenvalue_method(model, subcase)

        # op2 = self.op2
        # ---------------------------------------------------
        # title = ''
        # today = None
        # page_stamp = self.op2.make_stamp(title, today) # + '\n'

        node_gridtype = _get_node_gridtype(model, idtype=idtype)
        dof_map, ps = _get_dof_map(model)
        ngrid, ndof_per_grid, ndof_g = get_ndof(self.model, subcase)
        xyz_cid0 = model.grid.xyz_cid0()
        # -----------------------------------------------------------------------
        load_id, spc_id, mpc_id, suport_id, method_id = get_case_control(
            model.sol, subcase)
        dload_id, _ = subcase['DLOAD']
        freq_id, _ = subcase['FREQUENCY']
        # -----------------------------------------------------------------------

        GMN, mset = self.build_GMN(
            subcase, dof_map, ndof_g, xyz_cid0, fdtype=fdtype)
        self.GMN = GMN
        self.mset = mset
        mset_b = mset

        page_num, out = _run_modes(
            model, subcase,
            self.op2, f06_file,
            ngrid,
            ndof_per_grid,
            ndof_g,
            node_gridtype,
            dof_map,
            GMN=GMN,
            mset=mset,
            idtype=idtype,
            fdtype=fdtype,
            page_stamp=page_stamp, page_num=page_num,
            solver_dict=self.solver_dict,
            matrices_to_save={'Koa', 'Koo'},
        )
        #out['GMN'] = GMN
        #out['mset_bool'] = mset

        Koa = out.get('Koa')
        Koo = out.get('Koo')
        sset_b = out.get('sset_bool')
        oset_b = out.get('oset_bool')

        aset = out["aset"]
        Kaa = out["Kaa"]
        Maa = out["Maa"]
        eigenvalue = out["modes_eigenvalue"]
        omegas = np.sqrt(np.abs(eigenvalue))

        freqs = get_frequencies(model, subcase, omegas)
        nfreq = len(freqs)
        omegai = 2 * np.pi * freqs

        Fb = self.build_complex_Fb(model, subcase, ndof_g)
        Fg = xb_to_xg(model, Fb, ngrid, ndof_per_grid)
        Fa = reduce_Pg_to_Pa(
            Fg,
            ndof_g,
            sset_b,
            mset_b=mset_b,
            oset_b=oset_b,
            Gm=GMN, Koo=Koo, Koa=Koa)

        ndof_a = Kaa.shape[0]
        Caa = np.eye(ndof_a)
        np.fill_diagonal(Caa, 0.01)

        # nfreq = len(omegas)
        shape = (ndof_g, ndof_g, nfreq)
        tf_matrix = np.zeros(shape, dtype="complex128")  # full system 3x3 TF matrix
        tf_mask = np.ones((ndof_g, ndof_g, nfreq), dtype="bool")
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


def _save_eigenvectors(model: BDF,
                       subcase: Subcase,
                       phig: np.ndarray,
                       eigenvalue: np.ndarray,
                       mode_cycle: np.ndarray,
                       node_gridtype: np.ndarray,
                       op2: OP2,
                       f06_file: TextIO,
                       page_num: int,
                       sol_type: str,
                       page_stamp: str,
                       title: str, subtitle: str, label: str,
                       idtype: str='int32') -> int:
    isubcase = subcase.id
    nmode = len(eigenvalue)
    modes = np.arange(1, nmode + 1, dtype=idtype)
    eigenvalue = eigenvalue.astype("float32")
    (
        write_phi_f06,
        write_phi_op2,
        write_phi_pch,
        unused_options,
        phi_set,
    ) = get_f06_op2_pch_set(subcase, "DISPLACEMENT")
    write_eigenvector = any((write_phi_f06, write_phi_op2))

    if write_eigenvector:
        nnode = node_gridtype.shape[0]
        # print(xg_out.shape)
        node_gridtypei, phi, nnodei = slice_modal_set(
            node_gridtype, phig, nnode, nmode, phi_set)

        # phi:  (nmode, nnode*6)
        # data: (nmode, nnode, 6)
        data = phi.reshape((nmode, nnodei, 6)).astype("float32")
        table_name = "OUGV1"
        if sol_type == 'modes':
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
        else:
            assert sol_type == 'buckling', sol_type
            eigenvector_obj = RealEigenvectorArray.add_buckling_case(
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

    if write_phi_op2:
        op2.eigenvectors[isubcase] = eigenvector_obj

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
    return page_num

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


def get_rset(model: BDF,
             suport_id: int=0) -> set[tuple[int, int]]:
    """Creates the r-set from SUPORT/SUPORT1 cards.

    The r-set contains the reference DOFs used to define rigid-body motion
    in free-body (inertia relief) analysis.
    """
    rset_map: set[tuple[int, int]] = set()

    suport = model.suport
    if suport.n == 0:
        return rset_map

    suport_ids = []
    if 0 in suport.suport_id:
        suport_ids.append(0)
    if suport_id > 0:
        suport_ids.append(suport_id)
    suport = model.suport.slice_card_by_id(suport_id)

    for nid, comp in zip(suport.node_id, suport.component):
        comp_str = str(comp)
        for compi in comp_str:
            if compi != '0':
                rset_map.add((nid, int(compi)))
    return rset_map


def get_rset_bool(model: BDF,
                  dof_map: DOF_MAP,
                  ndof: int,
                  suport_id: int) -> np.ndarray:
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
    rset_map = get_rset(model, suport_id=suport_id)
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
    Kaa_dense = todense(Kaa)

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
        model: BDF,
        dof_map: DOF_MAP,
        fset: NDArrayNbool,
        suport_id: int=0,
        idtype: str = "int32") -> NDArrayNbool:
    """gets the residual structure dofs"""
    asetmap = get_aset(model)
    bsetmap = get_bset(model)
    csetmap = get_cset(model)
    rsetmap = get_rset(model, suport_id=suport_id)
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
                    "SUPORTi, CSETi, and BSETi entries.")

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
        set_map,
        dof_map: DOF_MAP,
        idtype: str = "int32",
        use_ints: bool = True) -> NDArrayNbool:
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


def xg_to_xb(model: BDF,
             xg: NDArrayNfloat,
             ngrid: int, ndof_per_grid: int,
             inplace: bool = True) -> NDArrayNfloat:
    assert isinstance(xg, np.ndarray)
    shape = xg.shape
    if xg.ndim == 1:
        ndof = len(xg)
        xg = xg.reshape(ndof, 1)
    else:
        assert xg.ndim == 2, xg.shape
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
            #xi = xg[i1:i2]   # 1d
            xi = xg[i1:i2, :]
            
            # TODO: verify the transform; I think it's right
            # (6,2) = (6,2) * (6,6,2)
            #xb[i1:i2] = xi @ T  # 1d
            xb[i1:i2, :] = xi @ T[:, :, np.newaxis]
    xb = xb.reshape(shape)
    return xb

def xb_to_xg(
    model: BDF,
    xb: NDArrayNfloat,
    ngrid: int, ndof_per_grid: int,
    inplace: bool = True) -> NDArrayNfloat:
    assert isinstance(xb, np.ndarray), type(xb)
    str(ngrid)

    xg = xb
    if not inplace:
        xg = copy.deepcopy(xg)

    nids = model._type_to_id_map["GRID"]
    for i, nid in enumerate(nids):
        node = model.nodes[nid]
        if node.cd:
            model.log.debug(f"node {nid} has a CD={node.cd}")
            cd_ref = node.cd_ref
            T = cd_ref.beta_n(n=2)
            i1 = i * ndof_per_grid
            i2 = (i + 1) * ndof_per_grid
            xi = xb[i1:i2]
            xg[i1:i2] = xi @ T.T  # TODO: verify the transform; I think it's right
    return xg


def write_oload_resultant(
    Fb: NDArrayNfloat,
    dof_map: DOF_MAP,
    isubcase: int,
    ngrid: int,
    ndof_per_grid: int,
    f06_file: TextIO,
    page_stamp: str,
    page_num: int,
    log: SimpleLogger,) -> int:
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
    log: SimpleLogger,) -> int:
    """Write SPCFORCE RESULTANT table to F06."""
    ndof = ngrid * ndof_per_grid
    #assert fspc.shape == (ndof,), (f'fspc.shape={fspc.shape}, ngrid={ngrid}, ngrid*6={ngrid*6}')
    fxyz_mxyz = fspc[:ngrid * ndof_per_grid].reshape(ngrid, ndof_per_grid)
    fxyz_mxyz_sum = fxyz_mxyz.sum(axis=0)
    spc_resultant = Resultant("SPCFORCE", fxyz_mxyz_sum, isubcase)
    log.info(f"SPCFORCE RESULTANT {fxyz_mxyz_sum}")
    spc_resultant.write_f06(f06_file, page_stamp, page_num)
    return page_num


def solve(Kaa: lil_matrix,
          Fa_solve: np.ndarray,
          aset: np.ndarray,
          log: SimpleLogger,
          idtype: str = "int32",):
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


def grid_point_weight(model: BDF, Mbb,
                      dof_map: DOF_MAP, ndof: int,
                      xyz_cid0: np.ndarray):
    xyz_cid0 = model.grid.xyz_cid0()
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
    cds = model.grid.cd
    cids = coord.coord_id
    # icds = coord.index(cds, assume_sorted=False, inverse=False)
    beta = coord.j
    for xyz0, cd in zip(xyz_cid0, cds):
        # coord = model.coord.s
        # print(f'nid={nid}')

        # TODO: what about otuput coordinate frames; specifically cylindrical and spherical frames?
        xi, yi, zi = xyz0 - dxyz
        Tr = np.array([
            [0, zi, -yi],
            [-zi, 0, xi],
            [yi, -xi, 0],
        ], dtype="float64",)
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


def _run_modes(
    model: BDF,
    subcase: Subcase,
    op2: OP2,
    f06_file,
    ngrid: int,
    ndof_per_grid: int,
    ndof: int,
    node_gridtype: np.ndarray,
    dof_map: dict,
    Kgg_override=None,
    Mgg_override=None,
    GMN=None,
    mset: np.ndarray | None = None,
    idtype: str = "int32",
    fdtype: str = "float64",
    title: str='', subtitle: str='', label: str='',
    page_stamp: str='',
    page_num: int=1,
    solver_dict=None,
    matrices_to_save: set[str] | None=None) -> tuple[int, dict[str, Any]]:
    if matrices_to_save is None:
        matrices_to_save = set([])
    assert solver_dict is not None, solver_dict

    neigenvalue, norm_str = get_real_eigenvalue_method(model, subcase)
    suport_id = subcase['SUPORT'][0] if 'SUPORT' in subcase else 0
    spc_id = subcase['SPC'][0] if 'SPC' in subcase else 0
    mpc_id = subcase['MPC'][0] if 'MPC' in subcase else 0

    log = model.log
    out = {}
    Kgg = get_Kgg(
        model, dof_map, ndof, ngrid, ndof_per_grid,
        idtype=idtype, fdtype=fdtype,
        Kgg_override=Kgg_override,
        solver_dict=solver_dict)
    out["Kgg"] = Kgg

    Mbb = build_Mbb(model, subcase, dof_map, ndof, fdtype=fdtype,
                    solver_dict=solver_dict)
    page_num = write_grid_point_weight(
        model, Mbb, dof_map, ndof,
        op2, f06_file,
        title=title, subtitle=subtitle, label=label,
        page_stamp=page_stamp, page_num=page_num)

    Mgg = Kbb_to_Kgg(model, Mbb, ngrid, ndof_per_grid)
    if 'Mgg' in matrices_to_save:
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
        #mset = None
        mset_b = None
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

    out['gset_bool'] = gset
    #out['gset_bool'] = gset_b
    out['sset'] = sset
    out['sset_bool'] = sset_b
    #out['mset'] = mset
    out['mset_bool'] = mset_b

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
    rset_b = get_rset_bool(model, dof_map, ndof, suport_id=suport_id)
    has_rset = np.any(rset_b) and GMN is None

    if has_rset:
        # r-set within the a-set
        a_indices = np.where(aset)[0]
        r_in_a = rset_b[a_indices]
        lset_local = np.where(~r_in_a)[0]
        rset_local = np.where(r_in_a)[0]
        nr = len(rset_local)
        nl = len(lset_local)
        log.warning(f"SUPORT r-set: {nr} DOFs, l-set: {nl} DOFs")

        Kaa_dense = todense(Kaa)
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
        phi_a = np.zeros((nmode, ndof_a), dtype=fdtype)
        for imode in range(nmode):
            phi_l = np.zeros(nl, dtype=fdtype)
            phi_l[ipositive] = xl_[:, imode]
            phi_a[imode, lset_local] = phi_l

        phi_n = np.full((nmode, ndof_work), 0.0, dtype=fdtype)
        phi_n[:, aset] = phi_a
        phi_n[:, sset] = xs
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
        phi_a = np.zeros((nmode, ndof_a), dtype=fdtype)
        for imode in range(nmode):
            phi_a[imode, ipositive] = xa_[:, imode]

        phi_n = np.full((nmode, ndof_work), 0.0, dtype=fdtype)
        phi_n[:, aset] = phi_a
        phi_n[:, sset] = xs

    phi_g = phi_n_to_g(phi_n, ndof, GMN=GMN, fdtype=fdtype)

    # phi_t = phi_g
    nnode_g = len(node_gridtype)
    phi, Mhh, Khh = apply_phi_normalization(
        Mgg, Kgg, eigenvalue, phi_g, nmode, nnode_g, norm_str)
    log.info(f"Mhh_diag: {np.diag(Mhh)}")
    log.info(f"Khh_diag: {np.diag(Khh)}")

    assert np.all(np.isfinite(phia))
    assert np.all(np.isfinite(phig))
    out["modes_eigenvalue"] = eigenvalue
    out["modes_norm"] = norm_str
    out["modes_phig"] = phig
    out["modes_Mhh"] = Mhh
    out["modes_Khh"] = Khh
    return page_num, out


def phi_n_to_g(phi_n: np.ndarray, ndof_g: int,
               GMN=None, fdtype: str='float64') -> np.ndarray:
    if GMN is None:
        return phi_n

    # ----------------------MPC recovery (n -> g) for modes----------------------
    nmode = phi_n.shape[1]

    # phig is currently in n-set (ndof_work = n_ndof), expand to g-set
    GMN_dense = GMN.toarray()
    phi_g = np.zeros((nmode, ndof_g), dtype=fdtype)
    for imode in range(nmode):
        phi_g[imode, :] = (GMN_dense @ phi_n[imode, :]).ravel()
    return phi_g


def _run_buckling(model: BDF,
                  subcase: Subcase,
                  ndof: int,
                  dof_map: DOF_MAP,
                  u_global: np.ndarray,
                  Kff,
                  free_dofs: np.ndarray) -> tuple[dict, int]:
    log = model.log
    neigenvalue, norm_str = get_real_eigenvalue_method(
        model, subcase)

    # Build geometric stiffness from preload stress state
    #KDgg = get_KDgg(
    #    model, dof_map, ndof,
    #    ngrid, ndof_per_grid,
    #    idtype=idtype, fdtype=fdtype,
    #    KDgg_override=self.KDgg_override)
    KDgg = build_KDgg(model, ndof, dof_map, u_global)
    KDff = KDgg.tocsc()[np.ix_(free_dofs, free_dofs)].toarray()

    log.debug("  geometric stiffness assembled")

    # Solve buckling eigenproblem: (K + lambda*KD)*x = 0
    # => K*x = -lambda*KD*x
    # Use scipy generalized eigenvalue: Kff @ x = lambda * (-KDff) @ x
    neg_KDff = -KDff

    # Use eigh for symmetric positive-definite B (standard buckling)
    # If -KD isn't positive definite, fall back to general eig
    try:
        eigenvalues_all, eigvecs_all = eigh(Kff, neg_KDff)
    except np.linalg.LinAlgError:
        eigenvalues_all, eigvecs_all = eig(Kff, neg_KDff)
        eigenvalues_all = eigenvalues_all.real
        eigvecs_all = eigvecs_all.real

    if 0:
        # Keep only positive eigenvalues
        # (physical buckling modes) sorted ascending
        pos_mask = eigenvalues_all > 0
        eigenvalues_pos = eigenvalues_all[pos_mask]
        eigvecs_pos = eigvecs_all[:, pos_mask]
    else:
        # user must limit the eigenvalues
        eigenvalues_pos = eigenvalues_all
        eigvecs_pos = eigvecs_all

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

    f06_file.write('\n                              M O D A L   E F F E C T I V E   M A S S\n')
    f06_file.write('\n')
    f06_file.write('  MODE    FREQUENCY        T1             T2             T3'
                   '             R1             R2             R3\n')

    for i in range(nmode):
        f06_file.write(f'  {i+1:4d}  {cycle[i]:12.6f}')
        for j in range(6):
            f06_file.write(f'  {eff_mass[i, j]:13.6E}')
        f06_file.write('\n')

    f06_file.write(
        '\n'
        '                      M O D A L   E F F E C T I V E   M A S S'
        '   F R A C T I O N\n\n')
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
    page_num: int = 1,) -> int:
    """Build grid point force balance and write to F06 via pyNastran object.

    For each grid, reports: applied load, SPC force, and the total
    (should be ~0 at equilibrium).
    """
    # Internal forces: F_int = K @ x (at each DOF)
    # f_internal = Kgg @ xg
    f_internal = dot2(Kgg, xg)

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
    X0: np.ndarray | None = None) -> tuple[np.ndarray, np.ndarray]:
    """
    Solve the generalized eigenproblem
      Kaa*x = lambda*Maa*x.

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
        eigenvalues, xa = scipy.linalg.eigh(Kaa2_dense, Maa2)
    elif use_lobpcg:
        X0_reduced = None
        if X0 is not None:
            X0_reduced = X0[is_modes, :]
        eigenvalues, xa_ = backend.lobpcg(Kaa2, k=neigenvalues, M=Maa2, X0=X0_reduced)
        xa = np.full((ndof, neigenvalues), np.nan, dtype=xa_.dtype)
        xa[no_modes, :] = 0
        xa[is_modes, :] = xa_
    else:
        try:
            eigenvalues, xa_ = backend.eigsh(
                Kaa2, k=neigenvalues, M=Maa2, which="SM", return_eigenvectors=True
            )
        except ArpackNoConvergence as e:
            raise FatalError(
                f'eigensolver did not converge: {e.args[0]}'
            ) from e
        xa = np.full((ndof, neigenvalues), np.nan, dtype=xa_.dtype)
        xa[no_modes, :] = 0
        xa[is_modes, :] = xa_
    return eigenvalues, xa


def _build_xg(model: BDF,
              dof_map: DOF_MAP, ndof: int,
              subcase: Subcase) -> NDArrayNfloat:
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
        #log.warning("no SPCs defined in case control")
        spc_set = np.array([], dtype="int32")
        sset = np.zeros(ndof, dtype="bool")
        return spc_set, sset, xspc

    log.warning('creating xg')
    spc_id, unused_options = subcase["SPC"]
 
    spcs = get_reduced_spcs(model, spc_id)
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


def get_reduced_spcs(model: BDF, spc_id: int) -> list:
    spcs = []
    if spc_id == 0:
        return spcs
    if model.spcadd.n > 0 and spc_id in model.spcadd:
        spc_ids = []
        spcadd = model.spcadd.slice_card_by_id(spc_id)
        #print(spcadd.get_stats())
        spc_ids = spcadd.spc_ids
    else:
        spc_ids = [spc_id]
        #raise RuntimeError('no spcs')
    #spcs = model.get_reduced_spcs(spc_id, consider_spcadd=True, stop_on_failure=True)

    spc_cards = [spc for spc in model.spc_cards if spc.n > 0]
    for spc_card in spc_cards:
        #print(spc_card.get_stats())
        for spc_idi in spc_ids:
            if spc_idi not in spc_card.spc_id:
                continue
            spci = spc_card.slice_card_by_id(spc_idi)
            spcs.append(spci)
    #assert len(spcs) > 0, spcs
    return spcs


def todense(matrix: np.ndarray):
    if hasattr(matrix, 'toarray'):
        dense_matrix = matrix.toarray()
    else:
        dense_matrix = np.asarray(matrix)
    return dense_matrix

def dot2(A, x):
    if hasattr(A, 'dot'):
        b = A.dot(x)
    else:
        b = A @ x
    return b


def _has_rigid_elements(model: BDF) -> bool:
    """Check if the model has any rigid elements (RBE2, RBE3, etc.)."""
    for attr in ("rbe2", "rbe3", "rbar", "rbar1", "rbe1", "rrod"):
        elem = getattr(model, attr, None)
        if elem is not None and elem.n > 0:
            return True
    return False


def _build_Gmn(model: BDF, mpc_id: int,
               dof_map: DOF_MAP, ndof: int,
               xyz_cid0: np.ndarray,
               fdtype: str = "float64",
               solver_dict=None,
               ) -> tuple[np.ndarray | None, np.ndarray | None]:
    """
    Supports:
     - MPC, MPCADD
     - RBE1, RBE2, RBE3
    """
    assert solver_dict is not None, solver_dict
    log = model.log

    # Collect all dependent DOFs from MPC cards and rigid elements
    dependents_list = []

    # --- Rigid elements (RBE2, RBE3, RBAR, RBAR1, RBE1, RROD) ---
    has_rigid = _has_rigid_elements(model)
    if mpc_id != 0 or has_rigid:
        log.info('creating Gmn')

    # --- MPC cards ---
    mpc = reduce_mpc_cards(model, mpc_id)
    for mpc_idi, (idim0, idim1) in zip(mpc.mpc_id, mpc.idim):
        nodes = mpc.node_id[idim0:idim1]
        components = mpc.components[idim0:idim1]
        nid_dep = nodes[0]
        comp_dep = int(components[0])
        idof_dep = dof_map[(nid_dep, comp_dep)]
        dependents_list.append(idof_dep)

    # --- Rigid elements (RBE2, RBE3, RBAR, RBAR1, RBE1, RROD) ---
    rigid_gmn_rows = None
    rigid_m_set = None

    if has_rigid:
        rigid_gmn_rows, rigid_m_set, _ = assemble_gmn(
            model, dof_map=dof_map, ndof=ndof, apply_cd=True)

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
    for mpc_idi, (idim0, idim1) in zip(mpc.mpc_id, mpc.idim):
        coefficients = mpc.coefficients[idim0:idim1]
        components = mpc.components[idim0:idim1]
        nodes = mpc.node_id[idim0:idim1]

        nid_dep = nodes[0]
        comp_dep = int(components[0])
        coeff_dep = coefficients[0]
        idof_dep = dof_map[(nid_dep, comp_dep)]

        for i in range(1, len(nodes)):
            nid_ind = nodes[i]
            comp_ind = int(components[i])
            coeff_ind = coefficients[i]
            if coeff_ind == 0.0:
                continue
            idof_ind = dof_map[(nid_ind, comp_ind)]
            n_col = g_to_n[idof_ind]
            if n_col < 0:
                log.warning(
                    f"MPC independent DOF ({nid_ind},{comp_ind}) is also "
                    "dependent in another constraint — skipping")
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

    write_mat(model, GMN, 'PRTGMN', solver_dict)

    GMN_csc = GMN.tocsc()
    n_mpc = len(mpc)
    n_rigid = len(rigid_m_set) if rigid_m_set is not None else 0
    log.info(f"MPC/rigid reduction: {len(mset)} dependent DOFs eliminated "
             f"(MPC={n_mpc}, rigid={n_rigid}, g={ndof} -> n={n_ndof})")
    return GMN_csc, mset


def reduce_mpc_cards(model: BDF,
                     mpc_id: int) -> MPC:
    """gets the MPCs; handles MPCADD"""
    if mpc_id == 0:
        return MPC(model)

    mpcadd = model.mpcadd
    mpc = model.mpc

    # handle MPCADD
    if mpc_id in mpcadd.mpc_id:
        mpcadd = mpcadd.slice_card_by_id(mpc_id)
        mpc_ids = mpcadd.ids
    else:
        impc = mpc.index(mpc_id)
        mpc_ids = [mpc_id]

    # handle MPCs
    mpc = mpc.slice_card_by_id(mpc_id, sort_ids=True)
    return mpc


def get_Kgg(model: BDF, dof_map: DOF_MAP,
            ndof: int, ngrid: int, ndof_per_grid: int,
            idtype: str='int32', fdtype: str='float64',
            Kgg_override=None,
            solver_dict=None):
    assert solver_dict is not None, solver_dict

    if Kgg_override is not None:
        # if issparse(Kgg_override):
        #     Kgg = Kgg_override
        # else:
        #     Kgg = csc_matrix(Kgg_override[:ndof, :ndof])
        if issparse(Kgg_override):
            Kgg = Kgg_override[:ndof, :ndof].tocsc()
        else:
            Kgg = csc_matrix(Kgg_override[:ndof, :ndof])
    else:
        Kgg = build_Kgg(
            model, dof_map, ndof, ngrid, ndof_per_grid,
            idtype=idtype, fdtype=fdtype)
        write_mat(model, Kgg, 'PRTKGG', solver_dict)
    assert Kgg is not None, Kgg
    return Kgg


def write_mat(model: BDF,
              Kgg,
              name: str,
              solver_dict: dict):
    assert isinstance(name, str), name
    assert isinstance(solver_dict, dict), solver_dict
    name2 = name[3:].upper()  # drop PRT

    prt_f06 = False
    prt_op4, prt_h5 = get_matprn(model, name)

    if Kgg is None:
        return
    
    if prt_f06 and 'f06' in solver_dict:
        solver_dict['f06'].write(name2, Kgg)
    if prt_op4 and 'op4' in solver_dict:
        solver_dict['op4'].write(name2, Kgg)
    if prt_h5 and 'h5' in solver_dict:
        solver_dict['h5'].write(name2, Kgg)


def get_matprn(model: BDF, name: str):
    prtmat = get_param(model, name, 0)
    # f06: 1
    # op2: 2
    # op4: 4
    # h5: 8
    prtmat = max(prtmat, 8+4+2+1)
    prt_h5 = prtmat // 8
    prtmat -= prt_h5 * 8

    prt_op4 = prtmat // 4
    prtmat -= prt_h5 * 4

    prt_op2 = prtmat // 2
    prtmat -= prt_op2 * 4

    prt_f06 = (prtmat == 1)
    return bool(prt_op4), bool(prt_h5)


def write_grid_point_weight(
        model: BDF, Mbb: Array | None,
        dof_map: DOF_MAP, ndof: int,
        op2: OP2,
        f06_file: TextIO,
        title: str='', subtitle: str='', label: str='',
        page_stamp: str='', page_num: int=1) -> int:
    if Mbb is not None:
        xyz_cid0 = model.grid.xyz_cid0()
        reference_point, MO = grid_point_weight(
            model, Mbb, dof_map, ndof, xyz_cid0=xyz_cid0)
        weight = make_grid_point_weight(
            reference_point,
            MO,
            approach_code=1,
            table_code=13,
            title=title,
            subtitle=subtitle,
            label=label,
            superelement_adaptivity_index="",)
        op2.grid_point_weight[label] = weight
        page_num = weight.write_f06(f06_file, page_stamp, page_num)
    return page_num


def recover_statics(model: BDF, op2: OP2,
                    f06_file: TextIO,
                    subcase: Subcase,
                    xg: np.ndarray,
                    Fg_oload: np.ndarray | None,
                    fspc: np.ndarray | None,
                    fmpc: np.ndarray | None,
                    node_gridtype: np.ndarray,
                    ngrid, ndof_per_grid,
                    title='', subtitle='', label='',
                    page_stamp: str='', page_num: int=1,
                    fdtype: str='float32') -> int:
    itime = 0
    ntimes = 1  # static
    #isubcase = subcase.id
    assert page_stamp is not None, page_stamp
    page_stamp % page_num
    page_num = _save_displacment(
        op2, f06_file,
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

    page_num = _save_applied_load(
        op2, f06_file,
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
        page_stamp=page_stamp)

    page_num = _save_spc_forces(
        op2, f06_file,
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
        page_stamp=page_stamp)

    page_num = _save_mpc_forces(
        op2, f06_file,
        subcase,
        itime,
        ntimes,
        node_gridtype,
        fmpc,
        ngrid,
        ndof_per_grid,
        title=title,
        subtitle=subtitle,
        label=label,
        fdtype=fdtype,
        page_num=page_num,
        page_stamp=page_stamp)

    return page_num


def get_cb_rset_dofs(model: BDF,
                     subcase: Subcase,
                     dof_map: DOF_MAP,
                     g_to_f: np.ndarray) -> list:
    """handles subcase"""
    if not len(model.suport):
        raise RuntimeError("SOL 31 requires SUPORT or SUPORT1 cards to define boundary DOFs")
    r_set_dofs = _get_rset_dofs(model, subcase, dof_map, g_to_f)
    return r_set_dofs


def get_arobc_set(model: BDF,
                  subcase: Subcase,
                  dof_map: DOF_MAP, g_to_f):
    # The a-set and o-set are created in the following ways:
    # 2. If ASETi or QSETi entries are present, then the a-set consists
    #    of all DOFs listed on ASETi entries and any entries listing its
    #    subsets, such as QSETi, SUPORTi, CSETi, and BSETi entries.  Any
    #    OMITi entries are redundant. The remaining f-set DOFs are placed
    #    in the o-set.

    r_set_dofs = _get_rset_dofs(model, subcase, dof_map, g_to_f)
    asetmap = get_aset(model)
    bsetmap = get_bset(model)
    csetmap = get_cset(model)
    osetmap = get_oset(model)
    qsetmap = get_qset(model)
    
    naset = len(asetmap)
    nbset = len(bsetmap)
    ncset = len(csetmap)
    noset = len(osetmap)
    nqset = len(qsetmap)
    nrset = len(r_set_dofs)
    if noset and not any(naset, nbset, ncset, nqset, nrset):
        # 1. If only OMITi entries are present, then the o-set consists
        #    of DOFs listed explicitly on OMITi entries. The remaining
        #    f-set DOFs are placed in the b-set, which is a subset of
        #    the a-set.
        del asetmap, bsetmap, csetmap, qsetmap, rsetmap
        # Koo
    elif not all(naset, nqset, noset) and any(nrset, nbset, ncset):
        # 3. If there are no ASETi, QSETi, or OMITi entries present but there
        #    are SUPORTi, BSETi, or CSETi entries present, then the entire
        #    f-set is placed in the a-set and the o-set is not created.
        del asetmap, qsetmap, osetmap
    elif any([naset, nqset]):
        # 2. If ASETi or QSETi entries are present, then the a-set consists
        #    of all DOFs listed on ASETi entries and any entries listing its
        #    subsets, such as QSETi, SUPORTi, CSETi, and BSETi entries.  Any
        #    OMITi entries are redundant. The remaining f-set DOFs are placed
        #    in the o-set.
        del osetmap
    else:
        raise RuntimeError('bad aset reduction')
      
    #aset_g = apply_dof_map_to_set(asetmap, dof_map, idtype=idtype)


def _get_rset_dofs(model: BDF,
                  subcase: Subcase,
                  dof_map: DOF_MAP,
                  g_to_f: np.ndarray) -> list:

    # handles subcase
    suport = model.suport
    
    suport_ids = []
    if 0 in suport.suport_id:
        suport_ids.append(0)

    if 'SUPORT' in subcase:
        suport_id, _ = subcase['SUPORT']
        suport_ids.append(suport_id)

    # loop over the suports
    r_set_dofs = []
    suport = suport.slice_card_by_id(suport_ids)
    for nid, comps in zip(suport.node_id, suport.component):
        for comp in str(comps):
            dof = int(comp)
            g_idx = dof_map[(nid, dof)]
            f_idx = g_to_f[g_idx]
            if f_idx >= 0:
                r_set_dofs.append((nid, dof))

    return r_set_dofs


def get_case_control(sol: int, subcase: Subcase):
    # case control checkout
    method_id = 0
    if sol == 101:
        load_id = subcase['LOAD'][0]
        spc_id = subcase['SPC'][0] if 'SPC' in subcase else 0
        mpc_id = subcase['MPC'][0] if 'MPC' in subcase else 0
        suport_id = subcase['SUPORT'][0] if 'SUPORT' in subcase else 0
    else:
        # modes/buckling/freq response
        load_id = subcase['LOAD'][0] if 'LOAD' in subcase else 0
        spc_id = subcase['SPC'][0] if 'SPC' in subcase else 0
        mpc_id = subcase['MPC'][0] if 'MPC' in subcase else 0
        method_id = subcase['METHOD'][0] if 'METHOD' in subcase else 0
        suport_id = subcase['SUPORT'][0] if 'SUPORT' in subcase else 0
    return load_id, spc_id, mpc_id, suport_id, method_id


def get_recover_results(subcase: Subcase,
                        results: list[str]) -> bool:
    recover_results = False
    for key in results:
        if key in subcase:
            value, options = subcase[key]
            if value != 'NONE':
                recover_results = True
                return recover_results
    return recover_results


def static_g_to_n(model: BDF,
                  Kgg, Mgg, GMN,
                  log: SimpleLogger, solver_dict):
    if GMN is not None:
        log.info('Applying GMN to K: n = g - m')
        # Transform to n-set: Knn = GMN^T @ Kgg @ GMN, Mnn = GMN^T @ Mgg @ GMN
        Knn = GMN.T @ Kgg @ GMN
        write_mat(model, Knn, 'PRTKNN', solver_dict)

        Mnn = GMN.T @ Mgg @ GMN if Mgg is not None else None
        write_mat(model, Mnn, 'PRTMNN', solver_dict)

        # AUTOSPC on N-set: detect singular DOFs after MPC reduction
        n_ndof_temp = Knn.shape[0]
        autospc_dofs_n = autospc_n_set(Knn, n_ndof_temp, log)
        if len(autospc_dofs_n) > 0:
            # Add these to the SPC set in the n-set coordinate system
            # They will be picked up during the SPC partition below
            autospc_n_dofs = autospc_dofs_n
        else:
            autospc_n_dofs = np.array([], dtype="int32")

        #un, ipos, ineg = solve(
        #    Knn, Fn, nset, log, idtype='int32')
    else:
        Knn = Kgg
        Mnn = Mgg
        autospc_n_dofs = np.array([], dtype="int32")
    return Knn, Mnn, autospc_n_dofs

def static_Pg_to_Pn(model: BDF,
                    Pg, GMN,
                    gset, mset,
                    log: SimpleLogger, solver_dict):
    if GMN is not None:
        # Transform force to n-set
        Pn = GMN.T @ Pg
        write_mat(model, Pn, 'PRTPN', solver_dict)
        # Transform enforced displacements to n-set
        # n-set indices: remove m-set DOFs from the indexing
        nset_indices = np.setdiff1d(gset, mset)
    else:
        Pn = Pg
        nset_indices = gset
    return Pn, nset_indices

def static_Knn_to_Kff(
    model: BDF,
    Knn, Mnn, Pn,
    xs,
    nset,
    log: SimpleLogger, solver_dict):
    """
    n = f + s
    
    {Fn} = {Ff*} = [Kff Kfs] {xf?}
           {Fs?}   [Kfs Kss] {xs*}
    
    {Ff} = [Kff] {xf} + [Kfs] {xs*}
    {Fs} = [Kfs] {xf} + [Kss] {xs*}
    """
    sset = np.isfinite(xs)
    xs_abs = np.abs(xs)
    if sset.sum():
        fset = nset - sset
        f_idx = np.ix_(fset)
        s_idx = np.ix_(sset)
        ff_idx = np.ix_(fset, fset)
        fs_idx = np.ix_(fset, sset)

        Kff = Knn.todense()[ff_idx].tocsc()
        if Mff is not None:
            Mff = Mnn.todense()[ff_idx].tocsc()
        if np.nansum(xs_abs) == 0.0:
            # {Ff} = [Kff] {xf} + [Kfs] {xs*}
            # {Ff} = [Kff] {xf}
            Ff = Fn[f_idx]
        else:
            # {Ff} = [Kff] {xf}
            # [Kff] {xf} = {Ff} - [Kfs] {xs}
            Ff0 = Fn[f_idx]
            xs0 = xs[s_idx]
            Ff = Ff0 - Kfs @ xs0
        xf, ipos, ineg = solve(
            Kff, Ff, fset, log, idtype='int32')
    else:
        Kff = Knn
        Mff = Mnn
        Pf = Pn
        uf = un
        fset = nset
    return Kff, Mff, Pf, uf, fset

def static_Kff_to_Kaa(
    model: BDF,
    Kaa, Mff, Pf,
    xs,
    nset,
    log: SimpleLogger, solver_dict):
    """
    f = a + o

    {Ff} = {Fa*} = [Kaa Kao] {xa?}
           {Fo*}   [Koa Koo] {xo?}
    
    {Fa*} = [Kaa] {xa} + [Kao] {xo?}
    {Fo*} = [Kao] {xa} + [Koo] {xo?}

    {Fa} = [Kaa] {xa} + [Kao] xo?}
    {xo} = [Koo]^-1 ({Fo} - [Kao] {xa})
    
    {Fa} = [Kaa] {xa} + [Kao] [Koo]^-1 ({Fo} - [Kao] {xa})
    {Fa} = [Kaa] {xa}
           + [Kao] [Koo]^-1 {Fo}
           - [Kao] [Koo]^-1 [Kao] {xa}

    {Fa} = [Kaa] {xa}
           + [Kao]   [Koo]^-1 {Fo}
           - [Koa]^T [Koo]^-1 [Kao] {xa}

    {Fa} = ([Kaa] - [Koa]^T [Koo]^-1 [Kao] ) + {xa}
           + [Kao] [Koo]^-1 {Fo}

    {Fa} - [Kao] [Koo]^-1 {Fo}
       = ([Kaa] - [Koa]^T [Koo]^-1 [Kao] ) + {xa}
           
    {Fa'}  = {Fa} - [Kao] [Koo]^-1 {Fo}
    {Kaa'} = [Kaa] - [Koa]^T [Koo]^-1 [Kao]
    {Fa'} = [Kaa'] + {xa}
    """
    pass
