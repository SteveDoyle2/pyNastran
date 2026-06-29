"""Tests for residual vector (RESVEC) computation.

Verification strategy:
1. Static completeness: u_modal with RESVEC must equal K^{-1}*F exactly.
2. Orthogonality: residual vectors must be M-orthogonal to all elastic modes.
3. Mass normalization: psi.T @ M @ psi = I.
4. Effective mass: sum of modal effective mass (including residuals) = total mass.
"""
import unittest
import numpy as np
from scipy.sparse import csc_matrix
from scipy.sparse.linalg import spsolve
from scipy.linalg import eigh


from pyNastran.dev.bdf_vectorized3.solver.utils_modes import ( 
    compute_mass_participation)
from pyNastran.dev.bdf_vectorized3.solver.residual_vectors import (
    compute_residual_vectors,
    build_inertia_load_vectors,)


def _make_spring_mass_chain(n: int = 10, k: float = 1000.0, m: float = 1.0):
    """Build a simple spring-mass chain (1 DOF per node, n nodes, fixed at node 0).

    Returns K, M at the free DOFs (nodes 1..n-1), so size is (n-1) x (n-1).
    """
    # Full stiffness (n x n tridiagonal)
    K_full = np.zeros((n, n))
    for i in range(n - 1):
        K_full[i, i] += k
        K_full[i + 1, i + 1] += k
        K_full[i, i + 1] -= k
        K_full[i + 1, i] -= k

    # Full mass (lumped diagonal)
    M_full = np.diag(np.full(n, m))

    # Apply SPC at node 0 (remove row/col 0)
    Kaa = K_full[1:, 1:]
    Maa = M_full[1:, 1:]
    return Kaa, Maa


def _solve_modes(Kaa, Maa, n_modes):
    """Solve the eigenvalue problem and return mass-normalized modes."""
    eigenvalues, phi = eigh(Kaa, Maa, subset_by_index=[0, n_modes - 1])
    # Mass-normalize
    for i in range(n_modes):
        gen_mass = phi[:, i] @ Maa @ phi[:, i]
        phi[:, i] /= np.sqrt(gen_mass)
    return eigenvalues, phi


class TestResidualVectors(unittest.TestCase):
    def test_static_completeness_dense(self):
        """With RESVEC, modal static response must exactly match K^{-1}*F."""
        Kaa, Maa = _make_spring_mass_chain(n=20, k=5000.0, m=2.0)
        n = Kaa.shape[0]  # 19 DOFs

        # Keep only 5 modes (significant truncation)
        n_modes = 5
        eigenvalues, phi = _solve_modes(Kaa, Maa, n_modes)

        # Apply a force at the tip (last DOF)
        F = np.zeros((n, 1))
        F[-1, 0] = 100.0

        # Exact static solution
        u_exact = np.linalg.solve(Kaa, F[:, 0])

        # Modal static (without RESVEC) — truncated, inaccurate
        phi_t_F = phi.T @ F[:, 0]  # (p,)
        u_modal_trunc = phi @ (phi_t_F / eigenvalues)

        # Compute residual vectors
        psi, omega2_pseudo = compute_residual_vectors(Kaa, Maa, phi, eigenvalues, F, tiny=1e-12)
        assert psi.shape[1] >= 1, f"Expected at least 1 residual vector, got {psi.shape[1]}"

        # Modal static WITH RESVEC
        phi_aug = np.hstack([phi, psi])
        eig_aug = np.concatenate([eigenvalues, omega2_pseudo])
        phi_aug_t_F = phi_aug.T @ F[:, 0]
        u_modal_aug = phi_aug @ (phi_aug_t_F / eig_aug)

        # Without RESVEC: significant error
        error_trunc = np.linalg.norm(u_modal_trunc - u_exact) / np.linalg.norm(u_exact)
        assert error_trunc > 0.01, f"Truncated modal should have >1% error, got {error_trunc:.6f}"

        # With RESVEC: machine precision
        error_aug = np.linalg.norm(u_modal_aug - u_exact) / np.linalg.norm(u_exact)
        assert error_aug < 1e-10, (
            f"Augmented modal should match static exactly, got error={error_aug:.2e}")


    def test_static_completeness_sparse(self):
        """Same as dense test but with sparse K."""
        Kaa_dense, Maa = _make_spring_mass_chain(n=30, k=8000.0, m=1.5)
        Kaa = csc_matrix(Kaa_dense)
        n = Kaa.shape[0]

        n_modes = 4
        eigenvalues, phi = _solve_modes(Kaa_dense, Maa, n_modes)

        # Force at mid-span
        F = np.zeros((n, 1))
        F[n // 2, 0] = 50.0

        u_exact = spsolve(Kaa, F[:, 0])

        psi, omega2_pseudo = compute_residual_vectors(Kaa, Maa, phi, eigenvalues, F, tiny=1e-12)

        phi_aug = np.hstack([phi, psi])
        eig_aug = np.concatenate([eigenvalues, omega2_pseudo])
        phi_aug_t_F = phi_aug.T @ F[:, 0]
        u_modal_aug = phi_aug @ (phi_aug_t_F / eig_aug)

        error = np.linalg.norm(u_modal_aug - u_exact) / np.linalg.norm(u_exact)
        assert error < 1e-10, f"Sparse RESVEC error={error:.2e}"


    def test_orthogonality(self):
        """Residual vectors must be M-orthogonal to all retained elastic modes."""
        Kaa, Maa = _make_spring_mass_chain(n=15, k=3000.0, m=1.0)
        n = Kaa.shape[0]

        n_modes = 6
        eigenvalues, phi = _solve_modes(Kaa, Maa, n_modes)

        # Multiple load vectors
        F = np.zeros((n, 3))
        F[0, 0] = 1.0
        F[n // 2, 1] = 1.0
        F[-1, 2] = 1.0

        psi, _ = compute_residual_vectors(Kaa, Maa, phi, eigenvalues, F, tiny=1e-12)

        if psi.shape[1] == 0:
            return

        # Check M-orthogonality: phi.T @ M @ psi should be zero
        cross = phi.T @ Maa @ psi
        assert np.allclose(cross, 0.0, atol=1e-12), (
            f"Cross M-orthogonality failed, max={np.max(np.abs(cross)):.2e}"
        )

        # Check self M-orthonormality of residual vectors
        psi_M_psi = psi.T @ Maa @ psi
        assert np.allclose(psi_M_psi, np.eye(psi.shape[1]), atol=1e-12), (
            f"Residual vectors not M-orthonormal, "
            f"max off-diag={np.max(np.abs(psi_M_psi - np.eye(psi.shape[1]))):.2e}"
        )


    def test_redundant_load_discarded(self):
        """A load already spanned by retained modes should produce no residual vector."""
        Kaa, Maa = _make_spring_mass_chain(n=10, k=1000.0, m=1.0)
        n = Kaa.shape[0]

        # Keep ALL modes (no truncation)
        n_modes = n
        eigenvalues, phi = _solve_modes(Kaa, Maa, n_modes)

        # Any load is fully spanned when all modes are retained
        F = np.zeros((n, 1))
        F[3, 0] = 42.0

        psi, _ = compute_residual_vectors(Kaa, Maa, phi, eigenvalues, F, tiny=1e-8)
        assert psi.shape[1] == 0, (
            f"With all modes retained, residual vector should be discarded. Got {psi.shape[1]}"
        )


    def test_effective_mass_completeness(self):
        """Sum of modal effective mass (with RESVEC) should equal total mass.

        For a 1D chain with unit mass per node, total mass = n-1.
        Modal effective mass in the chain direction = (phi_i.T @ M @ {1})^2.
        With RESVEC from ALL nodes as unit loads, the sum should equal total mass.

        The key insight: one inertia load (M*{1}) produces one residual vector,
        which captures one "direction" of missing mass. To capture ALL missing
        mass, we need residual vectors from enough independent loads to span
        the null space of the retained modal participation. Using unit forces
        at every DOF guarantees this.
        """
        n = 12
        Kaa, Maa = _make_spring_mass_chain(n=n, k=2000.0, m=1.0)
        ndof = Kaa.shape[0]  # n-1

        n_modes = 4
        eigenvalues, phi = _solve_modes(Kaa, Maa, n_modes)

        total_mass = np.sum(np.diag(Maa))  # = (n-1) * m

        # Direction vector (all translations in same direction for 1D)
        direction = np.ones(ndof)

        # Effective mass without RESVEC
        gamma = phi.T @ Maa @ direction
        eff_mass_trunc = np.sum(gamma**2)

        # Use unit force at every DOF as load vectors (guarantees completeness)
        F_all = np.eye(ndof)
        psi, omega2_pseudo = compute_residual_vectors(Kaa, Maa, phi, eigenvalues, F_all, tiny=1e-12)

        phi_aug = np.hstack([phi, psi])
        gamma_aug = phi_aug.T @ Maa @ direction
        eff_mass_aug = np.sum(gamma_aug**2)

        # Without RESVEC: less than total mass
        assert eff_mass_trunc < total_mass * 0.99, (
            f"Truncated effective mass should be < total. "
            f"Got {eff_mass_trunc:.4f} vs total {total_mass:.4f}"
        )

        # With RESVEC (all DOFs): equals total mass
        assert abs(eff_mass_aug - total_mass) < 1e-8, (
            f"Augmented effective mass should equal total mass. "
            f"Got {eff_mass_aug:.10f} vs {total_mass:.10f}, diff={abs(eff_mass_aug - total_mass):.2e}"
        )


    def test_multiple_loads_independent(self):
        """Multiple independent loads should produce multiple residual vectors."""
        Kaa, Maa = _make_spring_mass_chain(n=25, k=4000.0, m=1.0)
        n = Kaa.shape[0]

        n_modes = 3
        eigenvalues, phi = _solve_modes(Kaa, Maa, n_modes)

        # 4 independent point loads at different locations
        F = np.zeros((n, 4))
        F[2, 0] = 1.0
        F[8, 1] = 1.0
        F[14, 2] = 1.0
        F[20, 3] = 1.0

        psi, omega2_pseudo = compute_residual_vectors(Kaa, Maa, phi, eigenvalues, F, tiny=1e-12)

        # Should get multiple residual vectors (up to 4)
        assert psi.shape[1] >= 2, (
            f"Expected multiple residual vectors from 4 independent loads, got {psi.shape[1]}"
        )

        # Each residual vector should have positive pseudo-frequency
        assert np.all(omega2_pseudo > 0), "Pseudo-eigenvalues should be positive"

        # Pseudo-frequencies should be larger than highest retained mode
        assert np.all(omega2_pseudo > eigenvalues[-1]), (
            f"Pseudo-eigenvalues should exceed highest retained mode. "
            f"min pseudo={omega2_pseudo.min():.1f}, max retained={eigenvalues[-1]:.1f}"
        )


    def test_build_inertia_load_vectors(self):
        """Inertia load vectors should have correct structure."""
        n = 18  # 3 grids x 6 DOF
        Maa = np.diag(np.ones(n) * 2.0)

        F_inertia = build_inertia_load_vectors(Maa, ndof=6)

        assert F_inertia.shape == (n, 6)

        # TX direction: mass at DOFs 0, 6, 12
        expected_tx = np.zeros(n)
        expected_tx[0::6] = 2.0
        assert np.allclose(F_inertia[:, 0], expected_tx), "TX inertia load incorrect"

        # TY direction: mass at DOFs 1, 7, 13
        expected_ty = np.zeros(n)
        expected_ty[1::6] = 2.0
        assert np.allclose(F_inertia[:, 1], expected_ty), "TY inertia load incorrect"


    def test_drive_point_frf_at_dc(self):
        """Drive-point FRF at omega=0 with RESVEC should match K^{-1}(j,j).

        This is the classic RESVEC verification: the static compliance at the
        excitation point must be captured exactly.
        """
        Kaa, Maa = _make_spring_mass_chain(n=20, k=6000.0, m=1.0)
        n = Kaa.shape[0]
        j = n - 1  # drive point = tip

        n_modes = 5
        eigenvalues, phi = _solve_modes(Kaa, Maa, n_modes)

        # Static compliance at drive point
        K_inv = np.linalg.inv(Kaa)
        H_static_exact = K_inv[j, j]

        # Modal approximation of H(0) = sum phi_j^2 / omega_j^2
        H_trunc = np.sum(phi[j, :] ** 2 / eigenvalues)

        # With RESVEC
        F = np.zeros((n, 1))
        F[j, 0] = 1.0
        psi, omega2_pseudo = compute_residual_vectors(Kaa, Maa, phi, eigenvalues, F, tiny=1e-12)
        phi_aug = np.hstack([phi, psi])
        eig_aug = np.concatenate([eigenvalues, omega2_pseudo])
        H_aug = np.sum(phi_aug[j, :] ** 2 / eig_aug)

        # Truncated: error
        error_trunc = abs(H_trunc - H_static_exact) / H_static_exact
        assert error_trunc > 0.01, f"Expected truncation error, got {error_trunc:.6f}"

        # Augmented: exact
        error_aug = abs(H_aug - H_static_exact) / H_static_exact
        assert error_aug < 1e-10, f"Drive-point H(0) with RESVEC should be exact. Error={error_aug:.2e}"


    def test_mass_participation_factors(self):
        """Mass participation factors should sum to total mass per direction.

        For a 1D chain (single translational DOF per node), the sum of effective
        masses across ALL modes must equal the total translational mass.
        We also check that partial sums are monotonically increasing and that
        the function handles the 6-DOF-per-node layout correctly.
        """
        n = 8
        Kaa, Maa = _make_spring_mass_chain(n=n, k=5000.0, m=2.0)
        ndof = Kaa.shape[0]  # n-1 = 7

        # Solve ALL modes so effective mass sums exactly to total
        eigenvalues, phi_cols = eigh(Kaa, Maa)
        nmode = ndof

        # Mass-normalize
        for i in range(nmode):
            gm = phi_cols[:, i] @ Maa @ phi_cols[:, i]
            phi_cols[:, i] /= np.sqrt(gm)

        # Expand 1-DOF-per-node into 6-DOF-per-node layout for the function
        # Place chain DOF into the T1 (x-translation) slot
        nnode_g = ndof
        ndof_6 = nnode_g * 6
        Mgg = np.zeros((ndof_6, ndof_6))
        phig = np.zeros((nmode, ndof_6))

        for i in range(ndof):
            for j in range(ndof):
                Mgg[i * 6, j * 6] = Maa[i, j]
            for m in range(nmode):
                phig[m, i * 6] = phi_cols[i, m]

        mpf = compute_mass_participation(phig, Mgg, nnode_g, nmode)

        total_mass_t1 = mpf["total_mass"][0]
        sum_eff_mass_t1 = np.sum(mpf["effective_mass"][:, 0])

        # Sum of effective mass = total mass (all modes retained)
        assert abs(sum_eff_mass_t1 - total_mass_t1) < 1e-10, (
            f"sum eff mass={sum_eff_mass_t1:.10f}, total={total_mass_t1:.10f}"
        )

        # Cumulative ratio should reach 1.0 for T1
        assert abs(mpf["cumulative_ratio"][-1, 0] - 1.0) < 1e-10

        # Other directions (T2..R3) should have zero effective mass (no coupling)
        assert np.allclose(mpf["effective_mass"][:, 1:], 0.0, atol=1e-12)

        # Monotonically increasing cumulative ratio
        cum = mpf["cumulative_ratio"][:, 0]
        assert np.all(np.diff(cum) >= -1e-14)


    def test_mass_participation_partial_modes(self):
        """With fewer modes, cumulative ratio should be < 1.0."""
        n = 20
        Kaa, Maa = _make_spring_mass_chain(n=n, k=3000.0, m=1.5)
        ndof = Kaa.shape[0]

        n_modes = 5
        eigenvalues, phi_cols = _solve_modes(Kaa, Maa, n_modes)

        # Expand to 6-DOF layout
        nnode_g = ndof
        ndof_6 = nnode_g * 6
        Mgg = np.zeros((ndof_6, ndof_6))
        phig = np.zeros((n_modes, ndof_6))

        for i in range(ndof):
            for j in range(ndof):
                Mgg[i * 6, j * 6] = Maa[i, j]
            for m in range(n_modes):
                phig[m, i * 6] = phi_cols[i, m]

        mpf = compute_mass_participation(phig, Mgg, nnode_g, n_modes)

        # Partial modes: cumulative < 1.0
        assert mpf["cumulative_ratio"][-1, 0] < 1.0
        # But should capture significant mass (first mode of chain is dominant)
        assert mpf["cumulative_ratio"][-1, 0] > 0.5

        # Ratios should be non-negative
        assert np.all(mpf["effective_mass_ratio"] >= -1e-14)


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
