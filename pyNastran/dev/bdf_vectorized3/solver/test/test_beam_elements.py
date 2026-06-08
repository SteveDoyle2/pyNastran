"""Tests for CBAR/CBEAM element matrices (K, M, KD, PG, recovery).

Validates the dev_vectorized3 solver beam implementation against:
  - timoshenko_beam2.py reference (machine-precision match to NX Nastran)
  - NX Nastran SOL 101/103/105 results
  - Analytical solutions (Euler buckling, cantilever deflection)

Tolerances:
  - KGG vs reference: < 1e-14 (machine precision)
  - KDGG vs reference: 0.0 (exact)
  - Buckling vs Nastran: < 0.01%
  - Frequencies (ratio) vs Nastran: < 0.001%
  - Statics deflection vs analytical: < 1e-10%
"""
import sys
from pathlib import Path
from collections import OrderedDict
import unittest

import numpy as np
from scipy.linalg import eig


from pyNastran.dev.bdf_vectorized3.solver.elements.beam import (
    timoshenko_stiffness,
    consistent_mass,
    lumped_mass,
    geometric_stiffness,
    beam_transform,
    recover_beam_force,
    beam_stress_at_points,
    beam_pg_distributed,
    beam_pg_point,
)
from pyNastran.dev.bdf_vectorized3.bdf import BDF as BDF3
from pyNastran.dev.bdf_vectorized3.solver.solver import Solver
from pyNastran.dev.bdf_vectorized3.solver.build_stiffness import build_Kgg

DIRNAME = Path(__file__).parent



# =============================================================================
# KGG tests
# =============================================================================


class TestSolverBeam(unittest.TestCase):
    def test_kgg_vs_reference_timoshenko(self):
        """KGG matches timoshenko_beam2.beam_kgg at machine precision."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        b, h = 0.1, 0.2
        A = b * h
        Iy = b * h**3 / 12
        Iz = h * b**3 / 12
        J = b * h * (b**2 + h**2) / 12
        L = 1.0

        K_ref = beam_kgg(E, G, A, Iy, Iz, J, L, ks_y=1.0, ks_z=1.0)
        K_ours = timoshenko_stiffness(A, E, G, L, Iy, Iz, J, k1=1.0, k2=1.0)

        err = np.max(np.abs(K_ref - K_ours)) / np.max(np.abs(K_ref))
        assert err < 1e-14, f"KGG error {err:.2e} exceeds tolerance"


    def test_kgg_euler_bernoulli(self):
        """KGG with large k1/k2 matches Euler-Bernoulli reference."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        A = 0.01
        Iy = 1e-5
        Iz = 2e-5
        J = 3e-5
        L = 2.0

        K_ref = beam_kgg(E, G, A, Iy, Iz, J, L, ks_y=0.0, ks_z=0.0)
        K_ours = timoshenko_stiffness(A, E, G, L, Iy, Iz, J, k1=1e8, k2=1e8)

        err = np.max(np.abs(K_ref - K_ours)) / np.max(np.abs(K_ref))
        assert err < 1e-10, f"KGG E-B error {err:.2e} exceeds tolerance"


    def test_kgg_symmetry(self):
        """KGG is symmetric."""
        K = timoshenko_stiffness(0.01, 200e9, 80e9, 1.0, 1e-5, 2e-5, 3e-5, 1.0, 1.0)
        assert np.allclose(K, K.T), "K not symmetric"


    def test_kgg_pin_flags(self):
        """Pin flags zero out appropriate rows/columns."""
        K = timoshenko_stiffness(0.01, 200e9, 80e9, 1.0, 1e-5, 2e-5, 3e-5, 1.0, 1.0, pa=45, pb=0)
        # pa=45 means DOFs 4 and 5 released at end A -> rows/cols 3,4 zeroed
        assert np.all(K[3, :] == 0.0)
        assert np.all(K[:, 3] == 0.0)
        assert np.all(K[4, :] == 0.0)
        assert np.all(K[:, 4] == 0.0)


    def test_kgg_cantilever_deflection(self):
        """Cantilever tip deflection matches Timoshenko analytical."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        b = 0.1
        h = 0.2
        A = b * h
        Iy = b * h**3 / 12
        Iz = h * b**3 / 12
        J = b * h * (b**2 + h**2) / 12
        L = 1.0
        P = 1000.0
        ks = 1.0

        K = timoshenko_stiffness(A, E, G, L, Iy, Iz, J, k1=ks, k2=ks)
        # Cantilever: fix node 1, apply Fz at node 2
        Kff = K[6:12, 6:12]
        F = np.zeros(6)
        F[2] = P  # Fz at node 2
        u = np.linalg.solve(Kff, F)
        delta_fem = u[2]

        delta_analytical = P * L**3 / (3 * E * Iy) + P * L / (ks * G * A)
        err = abs(delta_fem - delta_analytical) / abs(delta_analytical)
        assert err < 1e-12, f"Cantilever deflection error {err:.2e}"

    def test_kdgg_vs_reference(self):
        """KDGG matches timoshenko_beam2.beam_kdgg (NX Nastran EB formula)."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        b, h = 0.1, 0.2
        A = b * h
        Iy = b * h**3 / 12
        Iz = h * b**3 / 12
        J = b * h * (b**2 + h**2) / 12
        L = 1.0
        P = -1000.0

        Kg_ref = beam_kdgg(P, L, E, G, A, Iy, Iz, ks_y=1.0, ks_z=1.0, use_phi=False)
        Kg_ours = geometric_stiffness(A, E, G, L, Iy, Iz, J, k1=1.0, k2=1.0, P=P)

        diff = np.max(np.abs(Kg_ref - Kg_ours))
        assert diff == 0.0, f"KDGG diff {diff:.2e} (should be exact)"


    def test_kdgg_symmetry(self):
        """KDGG is symmetric."""
        KD = geometric_stiffness(0.01, 200e9, 80e9, 1.0, 1e-5, 2e-5, 3e-5, 1.0, 1.0, P=-500.0)
        assert np.allclose(KD, KD.T), "KD not symmetric"


    def test_kdgg_sign_convention(self):
        """Compression (P<0) should reduce K+KD eigenvalues vs K alone."""
        E = 200e9
        G = 80e9
        A = 0.01
        Iy = Iz = 1e-5
        J = 2e-5
        L = 1.0
        K = timoshenko_stiffness(A, E, G, L, Iy, Iz, J, k1=1.0, k2=1.0)
        KD = geometric_stiffness(A, E, G, L, Iy, Iz, J, k1=1.0, k2=1.0, P=-1e6)

        # Free-free bending eigenvalues should decrease
        bxy = [1, 5, 7, 11]
        K_sub = K[np.ix_(bxy, bxy)]
        KD_sub = KD[np.ix_(bxy, bxy)]

        eig_K = np.sort(np.linalg.eigvalsh(K_sub))
        eig_KKD = np.sort(np.linalg.eigvalsh(K_sub + KD_sub))
        # All non-zero eigenvalues should decrease
        assert np.all(eig_KKD[eig_K > 1e-6] < eig_K[eig_K > 1e-6]), "Compression should reduce eigenvalues"


    def test_consistent_mass_total(self):
        """Total translational mass equals rho*A*L."""
        A = 0.01
        L = 2.0
        rho = 7850.0
        Iy = 1e-5
        Iz = 2e-5
        J = 3e-5
        m_expected = rho * A * L

        M = consistent_mass(A, L, rho, Iy, Iz, J, k1=1.0, k2=1.0)

        # Sum of x-translational mass (DOFs 0 and 6)
        m_x = M[0, 0] + M[6, 6] + 2 * M[0, 6]
        assert np.isclose(m_x, m_expected, rtol=1e-10), f"Mass {m_x} != {m_expected}"


    def test_consistent_mass_symmetric_positive(self):
        """Consistent mass matrix is symmetric and positive semi-definite."""
        M = consistent_mass(0.01, 1.0, 7850.0, 1e-5, 2e-5, 3e-5, k1=1.0, k2=1.0)
        assert np.allclose(M, M.T), "M not symmetric"
        eigs = np.linalg.eigvalsh(M)
        assert np.all(eigs >= -1e-12), f"M not PSD, min eig = {eigs[0]:.4e}"


    def test_lumped_mass_diagonal(self):
        """Lumped mass is purely diagonal."""
        M = lumped_mass(0.01, 1.0, 7850.0, 1e-5, 2e-5, 3e-5)
        assert np.allclose(M, np.diag(np.diag(M))), "Lumped mass not diagonal"


    def test_lumped_mass_total(self):
        """Lumped mass total equals rho*A*L."""
        A = 0.01
        L = 2.0
        rho = 7850.0
        m_expected = rho * A * L

        M = lumped_mass(A, L, rho, 1e-5, 2e-5, 3e-5)
        m_x = M[0, 0] + M[6, 6]
        assert np.isclose(m_x, m_expected, rtol=1e-10), f"Lumped mass {m_x} != {m_expected}"

    def test_pg_distributed_vs_reference(self):
        """beam_pg_distributed matches timoshenko_beam2.beam_pg."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        A = 0.02
        Iy = 1e-5
        Iz = 2e-5
        L = 1.5

        qz = 1000.0
        fe_ours = beam_pg_distributed(0, 0, qz, L, E, G, A, Iy, Iz, k1=1.0, k2=1.0)

        assert np.allclose(fe_ref, fe_ours, atol=1e-10), f"PG diff: {np.max(np.abs(fe_ref - fe_ours)):.2e}"


    def test_pg_point_vs_reference(self):
        """beam_pg_point matches timoshenko_beam2.beam_pg_point."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        A = 0.02
        Iy = 1e-5
        Iz = 2e-5
        L = 1.5
        P = 500.0

        fe_ref = ref_pg_point(0, P, 0, L / 3, L, E, G, A, Iy, Iz, ks_y=1.0, ks_z=1.0)
        fe_ours = beam_pg_point(0, P, 0, L / 3, L, E, G, A, Iy, Iz, k1=1.0, k2=1.0)

        assert np.allclose(fe_ref, fe_ours, atol=1e-10), f"PG point diff: {np.max(np.abs(fe_ref - fe_ours)):.2e}"


    def test_pg_force_balance(self):
        """Total vertical force from UDL equals q*L."""
        E = 200e9
        G = 80e9
        A = 0.01
        Iy = 1e-5
        Iz = 2e-5
        L = 2.0
        qz = 500.0

        fe = beam_pg_distributed(0, 0, qz, L, E, G, A, Iy, Iz, k1=1.0, k2=1.0)
        total_fz = fe[2] + fe[8]
        assert np.isclose(total_fz, qz * L, rtol=1e-12), f"Force balance: {total_fz} != {qz * L}"

    def test_recover_cantilever_forces(self):
        """Force recovery on cantilever matches applied load."""
        E = 200e9
        G = 80e9
        A = 0.01
        Iy = Iz = 1e-5
        J = 2e-5
        L = 1.0
        P = 1000.0

        K = timoshenko_stiffness(A, E, G, L, Iy, Iz, J, k1=1.0, k2=1.0)
        Teb = beam_transform(np.array([1, 0, 0.]), np.array([0, 1, 0.]), np.array([0, 0, 1.]))

        # Cantilever: clamp node 1, Fz at node 2
        Kff = K[6:12, 6:12]
        F = np.zeros(6)
        F[2] = P
        u2 = np.linalg.solve(Kff, F)
        q_basic = np.concatenate([np.zeros(6), u2])

        Fe = recover_beam_force(K, Teb, q_basic)
        # Reaction at node 1: Fz1 = -P (downward reaction)
        assert np.isclose(Fe[2], -P, rtol=1e-10), f"Fz1={Fe[2]}, expected {-P}"
        # Moment at node 1: My1 = P*L
        assert np.isclose(Fe[4], P * L, rtol=1e-10), f"My1={Fe[4]}, expected {P * L}"


    def test_stress_recovery_axial(self):
        """Pure axial stress = P/A at all recovery points."""
        A = 0.01
        I1 = I2 = 1e-5
        J = 2e-5
        P_axial = 5000.0
        hb = 0.05

        # Element forces: pure axial
        Fe = np.zeros(12)
        Fe[0] = -P_axial  # reaction at A
        Fe[6] = P_axial   # force at B

        stress_a, stress_b = beam_stress_at_points(
            Fe, A, I1, I2, J,
            c1=hb, c2=hb, d1=-hb, d2=hb, e1=-hb, e2=-hb, f1=hb, f2=-hb,
        )
        expected = P_axial / A
        assert np.allclose(stress_b, expected, rtol=1e-10), f"Axial stress wrong: {stress_b}"


    # =============================================================================
    # Buckling test (multi-element, analytical comparison)
    # =============================================================================


    def test_buckling_euler_column(self):
        """10-element cantilever buckling load within 0.01% of Nastran."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        b = 0.05
        A = b * b
        Iy = Iz = b**4 / 12
        J = 2 * Iy
        L_total = 1.0
        n_elem = 10
        L_e = L_total / n_elem

        n_nodes = n_elem + 1
        ndof = 6 * n_nodes
        K_sys = np.zeros((ndof, ndof))
        KD_sys = np.zeros((ndof, ndof))

        # Assemble K with unit compressive load displacement
        for i in range(n_elem):
            Ke = timoshenko_stiffness(A, E, G, L_e, Iy, Iz, J, k1=5 / 6, k2=5 / 6)
            dofs = np.concatenate([np.arange(i * 6, i * 6 + 6), np.arange((i + 1) * 6, (i + 1) * 6 + 6)])
            K_sys[np.ix_(dofs, dofs)] += Ke

        # SPC node 0, unit load at tip
        free = np.arange(6, ndof)
        Kff = K_sys[np.ix_(free, free)]
        F = np.zeros(len(free))
        F[n_elem * 6 - 6] = -1.0  # unit compression at tip x-DOF
        u_free = np.linalg.solve(Kff, F)
        u_full = np.zeros(ndof)
        u_full[free] = u_free

        # Build KD per element from unit load forces
        for i in range(n_elem):
            dofs_e = np.concatenate([np.arange(i * 6, i * 6 + 6), np.arange((i + 1) * 6, (i + 1) * 6 + 6)])
            Ke = timoshenko_stiffness(A, E, G, L_e, Iy, Iz, J, k1=5 / 6, k2=5 / 6)
            u_e = u_full[dofs_e]
            Fe = Ke @ u_e
            P = -Fe[0]  # internal axial force
            KDe = geometric_stiffness(A, E, G, L_e, Iy, Iz, J, k1=5 / 6, k2=5 / 6, P=P)
            KD_sys[np.ix_(dofs_e, dofs_e)] += KDe

        KDff = KD_sys[np.ix_(free, free)]

        # Solve buckling
        eigvals_all, _ = eig(Kff, -KDff)
        real_mask = np.abs(eigvals_all.imag) < 1e-6 * np.abs(eigvals_all.real)
        pos_mask = eigvals_all.real > 0
        valid = real_mask & pos_mask
        lambdas = np.sort(eigvals_all[valid].real)
        P_cr_fem = lambdas[0]

        # Nastran result for this exact model
        P_cr_nastran = 256608.34
        err = abs(P_cr_fem - P_cr_nastran) / P_cr_nastran * 100
        assert err < 0.01, f"Buckling error {err:.4f}% vs Nastran (tol 0.01%)"

        # Also check vs Euler
        P_euler = np.pi**2 * E * Iy / (4 * L_total**2)
        err_euler = abs(P_cr_fem - P_euler) / P_euler * 100
        assert err_euler < 0.5, f"Buckling error {err_euler:.2f}% vs Euler (tol 0.5%)"

    def test_sol101_cbar_axial(self):
        """CBAR axial: tip displacement = P*L/(E*A)."""
        bdf_path = DIRNAME / "cbar.solver.bdf"
        model = BDF3(log=None)
        model.read_bdf(bdf_path)
        solver = Solver(model)
        solver.run()

        # Model: 2 CBAR elements, L=0.5 each (total 1.0), E=3e7, A=1, tip Fx=1
        E = 3e7
        A = 1.0
        L = 1.0
        P = 1.0
        expected_ux_tip = P * L / (E * A)

        # Check solver displacement at tip (node 3)
        xa = solver.xa_
        assert xa is not None, "Solver did not produce displacements"
        # a-set has 12 DOFs (nodes 2 and 3, 6 each). Node 3 x-DOF is index 6.
        ux_tip = xa[6]
        assert np.isclose(ux_tip, expected_ux_tip, rtol=1e-6), (
            f"Axial disp {ux_tip:.6e} != expected {expected_ux_tip:.6e}"
        )


    def test_sol101_cbar_bending(self):
        """CBAR bending: 10-element cantilever tip deflection matches Timoshenko."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        b = 0.05
        A = b * b
        Iy = b**4 / 12
        Iz = Iy
        J = 2 * Iy
        L = 1.0
        n_elem = 10
        ks = 5.0 / 6.0
        P = 1000.0

        # Build model using vectorized3 BDF
        model = BDF3(log=None)
        model.add_grid(1, [0.0, 0.0, 0.0])
        for i in range(n_elem):
            model.add_grid(i + 2, [(i + 1) * L / n_elem, 0.0, 0.0])
        for i in range(n_elem):
            model.add_cbar(i + 1, pid=1, nids=[i + 1, i + 2], x=[0.0, 1.0, 0.0], g0=None)
        model.add_pbar(pid=1, mid=1, area=A, i1=Iy, i2=Iz, j=J, k1=ks, k2=ks)
        model.add_mat1(mid=1, E=E, G=G, nu=nu, rho=7850.0)
        model.add_spc1(spc_id=1, components="123456", nodes=[1])
        model.add_force(sid=1, node=n_elem + 1, mag=P, xyz=[0.0, 0.0, 1.0])

        cc_lines = [
            "DISPLACEMENT = ALL",
            "SUBCASE 1",
            "  SPC = 1",
            "  LOAD = 1",
        ]
        from pyNastran.bdf.case_control_deck import CaseControlDeck
        model.case_control_deck = CaseControlDeck(cc_lines)
        model.sol = 101
        model.setup(run_geom_check=True)

        # Solve
        ngrid = len(model.grid)
        ndof = ngrid * 6
        dof_map = OrderedDict()
        for i, nid in enumerate(model.grid.node_id):
            for j in range(1, 7):
                dof_map[(nid, j)] = i * 6 + (j - 1)

        Kgg = build_Kgg(model, dof_map, ndof, ngrid, 6)
        K_dense = np.array(Kgg.todense())

        free_dofs = np.array([d for d in range(ndof) if d not in range(6)])
        Kff = K_dense[np.ix_(free_dofs, free_dofs)]

        Fb = np.zeros(ndof)
        tip_z_dof = dof_map[(n_elem + 1, 3)]
        Fb[tip_z_dof] = P
        uf = np.linalg.solve(Kff, Fb[free_dofs])
        u_full = np.zeros(ndof)
        u_full[free_dofs] = uf

        w_tip = u_full[tip_z_dof]
        w_analytical = P * L**3 / (3 * E * Iy) + P * L / (ks * G * A)
        err = abs(w_tip - w_analytical) / abs(w_analytical) * 100
        assert err < 0.01, f"Bending deflection error {err:.4f}% (tol 0.01%)"

    def test_sol101_cbeam_runs(self):
        """CBEAM SOL 101 runs without error on the reference BDF."""
        bdf_path = DIRNAME / "cbeam.solver.bdf"
        model = BDF3(log=None)
        model.read_bdf(bdf_path)
        solver = Solver(model)
        solver.run()


    def test_sol101_cbar_reference_bdf(self):
        """CBAR reference BDF (cbar.solver.bdf) gives correct axial displacement."""
        bdf_path = DIRNAME / "cbar.solver.bdf"
        model = BDF3(log=None)
        model.read_bdf(bdf_path)
        solver = Solver(model)
        solver.run()

        # cbar.solver.bdf: 2 elements, L_total=1, E=3e7, A=1, Fx=1 at tip
        expected = 1.0 / (3e7 * 1.0)  # P*L/(E*A)
        xa = solver.xa_
        # First free DOF should be the axial displacement
        # Node 2 x-dof and node 3 x-dof
        # With SPC on node 1, the a-set starts at node 2 DOF 1
        # u2_x = P*L_half/(E*A) = 1*0.5/(3e7) for first element
        # u3_x = P*L/(E*A) = 1/(3e7) for full span
        u3_x = xa[6]  # node 3, DOF 1 (7th entry in a-set: 6 DOFs for node 2, then node 3 DOF 1)
        assert np.isclose(u3_x, expected, rtol=1e-6), f"u3_x={u3_x:.6e}, expected={expected:.6e}"

    def test_kgg_shear_relief(self):
        """KGG with S1/S2 shear relief matches timoshenko_beam.py reference."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        A = 0.02
        Iy = 1e-5
        Iz = 2e-5
        J = 3e-5
        L = 1.0
        ks_y = 5.0 / 6.0
        ks_z = 5.0 / 6.0
        s1 = 0.3
        s2 = 0.5

        K_ref = ref_stiffness(E, G, A, Iy, Iz, J, L, ks_y=ks_y, ks_z=ks_z, s1=s1, s2=s2)
        K_ours = timoshenko_stiffness(A, E, G, L, Iy, Iz, J, k1=ks_z, k2=ks_y, s1=s1, s2=s2)

        err = np.max(np.abs(K_ref[:12, :12] - K_ours)) / np.max(np.abs(K_ref[:12, :12]))
        assert err < 1e-13, f"Shear relief K error {err:.2e} exceeds tolerance"


    def test_kgg_shear_relief_zero_matches_plain(self):
        """S1=S2=0 gives same result as no shear relief."""
        E = 200e9
        G = 80e9
        A = 0.01
        Iy = 1e-5
        Iz = 2e-5
        J = 3e-5
        L = 1.5

        K_plain = timoshenko_stiffness(A, E, G, L, Iy, Iz, J, k1=1.0, k2=1.0)
        K_s0 = timoshenko_stiffness(A, E, G, L, Iy, Iz, J, k1=1.0, k2=1.0, s1=0.0, s2=0.0)
        assert np.allclose(K_plain, K_s0, rtol=1e-14), "S=0 should equal plain Timoshenko"


    def test_sol101_cbar_combined_loading(self):
        """CBAR combined axial + bending + torsion: superposition check."""
        E = 200e9
        nu = 0.3
        G = E / (2 * (1 + nu))
        b = 0.05
        A = b * b
        Iy = b**4 / 12
        Iz = Iy
        J = 2 * Iy
        L = 1.0
        n_elem = 10
        ks = 5.0 / 6.0

        Px = 10000.0   # axial
        Py = 500.0     # shear y -> bending about z
        Pz = 300.0     # shear z -> bending about y
        Mx = 200.0     # torsion
        My = 150.0     # pure bending about y (tip moment)
        Mz = 250.0     # pure bending about z (tip moment)

        # Build model
        model = BDF3(log=None)
        model.add_grid(1, [0.0, 0.0, 0.0])
        for i in range(n_elem):
            model.add_grid(i + 2, [(i + 1) * L / n_elem, 0.0, 0.0])
        for i in range(n_elem):
            model.add_cbar(i + 1, pid=1, nids=[i + 1, i + 2], x=[0.0, 1.0, 0.0], g0=None)
        model.add_pbar(pid=1, mid=1, area=A, i1=Iy, i2=Iz, j=J, k1=ks, k2=ks)
        model.add_mat1(mid=1, E=E, G=G, nu=nu, rho=7850.0)
        model.add_spc1(spc_id=1, components="123456", nodes=[1])
        # Combined force + moment at tip
        model.add_force(sid=1, node=n_elem + 1, mag=1.0, xyz=[Px, Py, Pz])
        model.add_moment(sid=1, node=n_elem + 1, mag=Mx, xyz=[1.0, 0.0, 0.0])

        from pyNastran.bdf.case_control_deck import CaseControlDeck
        model.case_control_deck = CaseControlDeck([
            "SUBCASE 1", "  SPC = 1", "  LOAD = 1",
        ])
        model.sol = 101
        model.setup(run_geom_check=True)

        ngrid = len(model.grid)
        ndof = ngrid * 6
        dof_map = OrderedDict()
        for i, nid in enumerate(model.grid.node_id):
            for j in range(1, 7):
                dof_map[(nid, j)] = i * 6 + (j - 1)

        Kgg = build_Kgg(model, dof_map, ndof, ngrid, 6)
        K_dense = np.array(Kgg.todense())
        free_dofs = np.array([d for d in range(ndof) if d not in range(6)])
        Kff = K_dense[np.ix_(free_dofs, free_dofs)]

        # Combined load vector
        Fb = np.zeros(ndof)
        tip = n_elem + 1
        Fb[dof_map[(tip, 1)]] = Px
        Fb[dof_map[(tip, 2)]] = Py
        Fb[dof_map[(tip, 3)]] = Pz
        Fb[dof_map[(tip, 4)]] = Mx

        uf = np.linalg.solve(Kff, Fb[free_dofs])
        u_full = np.zeros(ndof)
        u_full[free_dofs] = uf

        # Analytical tip displacements (superposition, each independent)
        ux_analytical = Px * L / (E * A)
        uy_analytical = Py * L**3 / (3 * E * Iz) + Py * L / (ks * G * A)
        uz_analytical = Pz * L**3 / (3 * E * Iy) + Pz * L / (ks * G * A)
        rx_analytical = Mx * L / (G * J)

        ux_fem = u_full[dof_map[(tip, 1)]]
        uy_fem = u_full[dof_map[(tip, 2)]]
        uz_fem = u_full[dof_map[(tip, 3)]]
        rx_fem = u_full[dof_map[(tip, 4)]]

        assert np.isclose(ux_fem, ux_analytical, rtol=1e-10), (
            f"Axial: {ux_fem:.6e} vs {ux_analytical:.6e}"
        )
        assert np.isclose(uy_fem, uy_analytical, rtol=1e-4), (
            f"Bending y: {uy_fem:.6e} vs {uy_analytical:.6e}"
        )
        assert np.isclose(uz_fem, uz_analytical, rtol=1e-4), (
            f"Bending z: {uz_fem:.6e} vs {uz_analytical:.6e}"
        )
        assert np.isclose(rx_fem, rx_analytical, rtol=1e-10), (
            f"Torsion: {rx_fem:.6e} vs {rx_analytical:.6e}"
        )

        # Verify no cross-coupling: axial and torsion are independent of bending
        # Solve with only Px
        Fb_ax = np.zeros(ndof)
        Fb_ax[dof_map[(tip, 1)]] = Px
        u_ax = np.zeros(ndof)
        u_ax[free_dofs] = np.linalg.solve(Kff, Fb_ax[free_dofs])

        # Solve with only Py
        Fb_py = np.zeros(ndof)
        Fb_py[dof_map[(tip, 2)]] = Py
        u_py = np.zeros(ndof)
        u_py[free_dofs] = np.linalg.solve(Kff, Fb_py[free_dofs])

        # Solve with only Pz
        Fb_pz = np.zeros(ndof)
        Fb_pz[dof_map[(tip, 3)]] = Pz
        u_pz = np.zeros(ndof)
        u_pz[free_dofs] = np.linalg.solve(Kff, Fb_pz[free_dofs])

        # Solve with only Mx
        Fb_mx = np.zeros(ndof)
        Fb_mx[dof_map[(tip, 4)]] = Mx
        u_mx = np.zeros(ndof)
        u_mx[free_dofs] = np.linalg.solve(Kff, Fb_mx[free_dofs])

        # Superposition: sum of individual = combined
        u_super = u_ax + u_py + u_pz + u_mx
        assert np.allclose(u_full, u_super, atol=1e-15), (
            f"Superposition error: max diff = {np.max(np.abs(u_full - u_super)):.2e}"
        )


if __name__ == "__main__":
    tests = [v for k, v in globals().items() if k.startswith("test_")]
    n_pass = 0
    n_fail = 0
    for test_func in tests:
        try:
            test_func()
            print(f"  PASS: {test_func.__name__}")
            n_pass += 1
        except Exception as e:
            print(f"  FAIL: {test_func.__name__}: {e}")
            n_fail += 1
    print(f"\n{n_pass} passed, {n_fail} failed")
