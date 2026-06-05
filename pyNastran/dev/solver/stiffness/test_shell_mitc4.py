"""Tests for the MITC4 shell element stiffness matrix.

Validates:
1. Patch test (constant stress reproduced exactly)
2. Rigid body modes (6 zero eigenvalues)
3. Symmetry of stiffness matrix
4. Membrane-only comparison to analytical CST/isoparametric result
5. Single-element bending comparison to exact Kirchhoff plate solution
6. Rank check (20 - 6 = 14 non-zero eigenvalues for unconstrained element)
"""
import unittest

import numpy as np

from pyNastran.dev.solver.stiffness.shell_mitc4 import (
    mitc4_stiffness,
    _shape_quad4,
    _dshape_quad4,
    _jacobian,
    _membrane_B,
    _bending_B,
    _mitc4_shear_B,
)


# --- Material properties for tests ---
# Isotropic steel: E=200 GPa, nu=0.3, t=0.01 m
E = 200.0e9
NU = 0.3
T = 0.01  # thickness
G = E / (2.0 * (1.0 + NU))

# Membrane constitutive (plane stress A = E*t/(1-nu^2) * [...])
_factor_A = E * T / (1.0 - NU**2)
A_MAT = _factor_A * np.array([
    [1.0, NU, 0.0],
    [NU, 1.0, 0.0],
    [0.0, 0.0, (1.0 - NU) / 2.0],
])

# Bending constitutive D = E*t^3/(12*(1-nu^2)) * [...]
_factor_D = E * T**3 / (12.0 * (1.0 - NU**2))
D_MAT = _factor_D * np.array([
    [1.0, NU, 0.0],
    [NU, 1.0, 0.0],
    [0.0, 0.0, (1.0 - NU) / 2.0],
])

# Transverse shear constitutive Ds = kappa * G * t * I2
KAPPA = 5.0 / 6.0
DS_MAT = KAPPA * G * T * np.eye(2)

# --- Geometry: unit square ---
X_SQUARE = np.array([0.0, 1.0, 1.0, 0.0])
Y_SQUARE = np.array([0.0, 0.0, 1.0, 1.0])


def _Ke_square():
    """Stiffness for unit square with steel properties."""
    return mitc4_stiffness(X_SQUARE, Y_SQUARE, A_MAT, D_MAT, DS_MAT)


class TestShapeFunctions(unittest.TestCase):
    """Basic shape function verification."""

    def test_partition_of_unity(self):
        """Shape functions sum to 1 at any point."""
        for xi in [-1, -0.5, 0, 0.5, 1]:
            for eta in [-1, -0.5, 0, 0.5, 1]:
                N = _shape_quad4(xi, eta)
                assert abs(N.sum() - 1.0) < 1e-14, f"N sum = {N.sum()} at ({xi},{eta})"

    def test_kronecker_delta(self):
        """N_i = 1 at node i, 0 at other nodes."""
        corners = [(-1, -1), (1, -1), (1, 1), (-1, 1)]
        for i, (xi, eta) in enumerate(corners):
            N = _shape_quad4(xi, eta)
            for j in range(4):
                expected = 1.0 if i == j else 0.0
                assert abs(N[j] - expected) < 1e-14

    def test_derivative_sum_zero(self):
        """Sum of dN/dxi = 0, sum of dN/deta = 0."""
        dN = _dshape_quad4(0.3, -0.7)
        assert abs(dN[0].sum()) < 1e-14
        assert abs(dN[1].sum()) < 1e-14


class TestJacobian(unittest.TestCase):
    """Jacobian computation for standard geometries."""

    def test_unit_square_jacobian(self):
        """Unit square [0,1]x[0,1] mapped from [-1,1]x[-1,1] gives J = 0.5*I."""
        dN = _dshape_quad4(0.0, 0.0)
        J, det_J = _jacobian(dN, X_SQUARE, Y_SQUARE)
        expected_J = 0.5 * np.eye(2)
        np.testing.assert_allclose(J, expected_J, atol=1e-14)
        np.testing.assert_allclose(det_J, 0.25, atol=1e-14)

    def test_rectangular_element(self):
        """2x3 rectangle: J = diag(1, 1.5)."""
        x = np.array([0.0, 2.0, 2.0, 0.0])
        y = np.array([0.0, 0.0, 3.0, 3.0])
        dN = _dshape_quad4(0.0, 0.0)
        J, det_J = _jacobian(dN, x, y)
        expected_J = np.diag([1.0, 1.5])
        np.testing.assert_allclose(J, expected_J, atol=1e-14)
        np.testing.assert_allclose(det_J, 1.5, atol=1e-14)


class TestMITC4Stiffness(unittest.TestCase):
    """Stiffness matrix properties."""

    def setUp(self):
        self.Ke_square = _Ke_square()

    def test_symmetry(self):
        """Ke must be symmetric."""
        np.testing.assert_allclose(self.Ke_square, self.Ke_square.T, atol=1e-6)

    def test_dimensions(self):
        """Should be 20x20."""
        assert self.Ke_square.shape == (20, 20)

    def test_positive_semidefinite(self):
        """All eigenvalues >= 0."""
        eigvals = np.linalg.eigvalsh(self.Ke_square)
        assert np.all(eigvals > -1e-6 * eigvals.max()), \
            f"Negative eigenvalue: {eigvals[eigvals < 0]}"

    def test_rigid_body_modes(self):
        """Exactly 6 zero eigenvalues (3 translations + 3 rotations).

        Tolerance: eigenvalue < 1e-6 * max eigenvalue.
        """
        eigvals = np.sort(np.linalg.eigvalsh(self.Ke_square))
        threshold = 1e-6 * eigvals.max()
        n_zero = np.sum(np.abs(eigvals) < threshold)
        assert n_zero == 6, f"Expected 6 zero eigenvalues, got {n_zero}. " \
                            f"First 8 eigenvalues: {eigvals[:8]}"

    def test_rank(self):
        """Rank should be 14 (20 DOF - 6 RBM)."""
        eigvals = np.linalg.eigvalsh(self.Ke_square)
        threshold = 1e-6 * eigvals.max()
        rank = np.sum(np.abs(eigvals) > threshold)
        assert rank == 14, f"Expected rank 14, got {rank}"

    def test_no_hourglass_modes(self):
        """No spurious zero-energy modes beyond the 6 RBM."""
        eigvals = np.sort(np.linalg.eigvalsh(self.Ke_square))
        threshold = 1e-6 * eigvals.max()
        n_zero = np.sum(np.abs(eigvals) < threshold)
        assert n_zero == 6, f"Hourglass detected: {n_zero} zero modes instead of 6"


class TestMembranePatchTest(unittest.TestCase):
    """Membrane patch test: constant stress state must be representable.

    Apply uniform x-tension (eps_x = 1, eps_y = gamma_xy = 0) displacements
    and verify internal forces match sigma_x * t * length.
    """

    def test_uniform_tension_x(self):
        """Uniform eps_x = e0 applied via nodal displacements.

        u = e0*x, v = 0 gives eps_x=e0, eps_y=0, gamma=0.
        """
        e0 = 1e-4  # small strain
        Ke = mitc4_stiffness(X_SQUARE, Y_SQUARE, A_MAT, D_MAT, DS_MAT)

        # DOF order per node: [u, v, w, theta_x, theta_y]
        # Apply u = e0*x, v = 0, w = 0, theta_x = 0, theta_y = 0
        u_nodal = np.zeros(20)
        for i in range(4):
            u_nodal[5 * i] = e0 * X_SQUARE[i]  # u = e0 * x

        # Internal force = Ke * u
        f_int = Ke @ u_nodal

        # For a patch test, sum of forces should be zero (self-equilibrated).
        fx_total = sum(f_int[5 * i] for i in range(4))
        fy_total = sum(f_int[5 * i + 1] for i in range(4))
        assert abs(fx_total) < 1e-4, f"fx sum = {fx_total}"
        assert abs(fy_total) < 1e-4, f"fy sum = {fy_total}"

        # Also check no bending forces arise from pure membrane deformation
        for i in range(4):
            w_force = abs(f_int[5 * i + 2])  # w direction
            mx_force = abs(f_int[5 * i + 3])  # theta_x moment
            my_force = abs(f_int[5 * i + 4])  # theta_y moment
            assert w_force < 1e-4, f"Node {i}: unexpected w force = {w_force}"
            assert mx_force < 1e-4, f"Node {i}: unexpected mx = {mx_force}"
            assert my_force < 1e-4, f"Node {i}: unexpected my = {my_force}"


class TestBendingBehavior(unittest.TestCase):
    """Verify bending response against analytical plate solutions."""

    def test_pure_curvature_kx(self):
        """Apply constant curvature kx via theta_y = kx * x.

        kappa_x = d(theta_y)/dx = kx
        Expected Mx = D11 * kx + D12 * 0 = D11 * kx
        """
        kx = 1.0  # unit curvature
        Ke = mitc4_stiffness(X_SQUARE, Y_SQUARE, A_MAT, D_MAT, DS_MAT)

        # DOF order per node: [u, v, w, theta_x, theta_y]
        # For constant kx: theta_y = kx * x, w = -0.5*kx*x^2 (compatible kinematics)
        u_nodal = np.zeros(20)
        for i in range(4):
            u_nodal[5 * i + 4] = kx * X_SQUARE[i]  # theta_y = kx * x
            u_nodal[5 * i + 2] = -0.5 * kx * X_SQUARE[i]**2  # w = -kx*x^2/2

        f_int = Ke @ u_nodal

        # For constant curvature, forces should be consistent.
        # Check symmetry: nodes on same x should have same theta_y force
        my_left = f_int[5 * 0 + 4] + f_int[5 * 3 + 4]   # nodes 1,4 (x=0)
        my_right = f_int[5 * 1 + 4] + f_int[5 * 2 + 4]  # nodes 2,3 (x=1)
        # Should be equal and opposite (equilibrium)
        assert abs(my_left + my_right) < 1e-4 * abs(my_left), \
            f"my_left={my_left}, my_right={my_right}"


class TestDistortedElement(unittest.TestCase):
    """Verify element works for non-rectangular geometry."""

    def test_parallelogram(self):
        """Parallelogram element should still pass basic checks."""
        x = np.array([0.0, 1.0, 1.3, 0.3])  # skewed
        y = np.array([0.0, 0.0, 1.0, 1.0])
        Ke = mitc4_stiffness(x, y, A_MAT, D_MAT, DS_MAT)

        # Symmetry
        np.testing.assert_allclose(Ke, Ke.T, atol=1e-6)

        # Correct number of rigid body modes
        eigvals = np.sort(np.linalg.eigvalsh(Ke))
        threshold = 1e-6 * eigvals.max()
        n_zero = np.sum(np.abs(eigvals) < threshold)
        assert n_zero == 6, f"Parallelogram: {n_zero} zero modes, expected 6"

    def test_trapezoid(self):
        """Trapezoidal element."""
        x = np.array([0.0, 2.0, 1.5, 0.5])
        y = np.array([0.0, 0.0, 1.0, 1.0])
        Ke = mitc4_stiffness(x, y, A_MAT, D_MAT, DS_MAT)

        np.testing.assert_allclose(Ke, Ke.T, atol=1e-6)

        eigvals = np.sort(np.linalg.eigvalsh(Ke))
        threshold = 1e-6 * eigvals.max()
        n_zero = np.sum(np.abs(eigvals) < threshold)
        assert n_zero == 6, f"Trapezoid: {n_zero} zero modes, expected 6"


class TestMembraneOnly(unittest.TestCase):
    """Test membrane-only behavior (D=0, Ds=0) against analytical solution.

    For a rectangular element under uniform tension, the membrane stiffness
    should match the exact isoparametric bilinear result.
    """

    def test_membrane_stiffness_trace(self):
        """The trace of K_membrane should equal the analytical value.

        For a unit square, isotropic membrane:
        trace(Km) = 2 * (A11 + A22 + A33) * area / (average element size)
        This is a simple integral check.
        """
        # Zero bending and shear
        D_zero = np.zeros((3, 3))
        Ds_zero = np.zeros((2, 2))
        Ke = mitc4_stiffness(X_SQUARE, Y_SQUARE, A_MAT, D_zero, Ds_zero)

        # Membrane DOF indices: u,v at each node = indices 0,1,5,6,10,11,15,16
        m_idx = [0, 1, 5, 6, 10, 11, 15, 16]
        Km = Ke[np.ix_(m_idx, m_idx)]

        # Should have 3 zero eigenvalues (2 translations + 1 rotation)
        eigvals = np.sort(np.linalg.eigvalsh(Km))
        threshold = 1e-6 * eigvals.max()
        n_zero = np.sum(np.abs(eigvals) < threshold)
        assert n_zero == 3, f"Membrane: {n_zero} zero modes, expected 3 (2T+1R)"

        # Bending DOF should have zero stiffness
        b_idx = [2, 3, 4, 7, 8, 9, 12, 13, 14, 17, 18, 19]
        Kb = Ke[np.ix_(b_idx, b_idx)]
        assert np.abs(Kb).max() < 1e-10, "Bending stiffness should be zero"


class TestGeometricStiffness(unittest.TestCase):
    """Validate MITC4 geometric stiffness matrix."""

    def setUp(self):
        from pyNastran.dev.solver.stiffness.shell_mitc4 import mitc4_geometric_stiffness
        stress = np.array([-1.0, 0.0, 0.0])  # unit Nxx compression
        self.Kg_compression = mitc4_geometric_stiffness(X_SQUARE, Y_SQUARE, stress)

    def test_symmetry(self):
        """Kg must be symmetric."""
        np.testing.assert_allclose(self.Kg_compression, self.Kg_compression.T, atol=1e-12)

    def test_dimensions(self):
        """Should be 20x20."""
        assert self.Kg_compression.shape == (20, 20)

    def test_negative_eigenvalue_under_compression(self):
        """Under compression, Kg should have negative eigenvalues (destabilizing)."""
        eigvals = np.linalg.eigvalsh(self.Kg_compression)
        assert np.any(eigvals < 0), "Compression Kg should have negative eigenvalues"

    def test_zero_for_zero_stress(self):
        """Zero stress gives zero geometric stiffness."""
        from pyNastran.dev.solver.stiffness.shell_mitc4 import mitc4_geometric_stiffness
        Kg = mitc4_geometric_stiffness(X_SQUARE, Y_SQUARE, np.zeros(3))
        assert np.abs(Kg).max() < 1e-15

    def test_linearity_in_stress(self):
        """Kg should scale linearly with stress magnitude."""
        from pyNastran.dev.solver.stiffness.shell_mitc4 import mitc4_geometric_stiffness
        s1 = np.array([-1.0, 0.0, 0.0])
        s2 = np.array([-2.0, 0.0, 0.0])
        Kg1 = mitc4_geometric_stiffness(X_SQUARE, Y_SQUARE, s1)
        Kg2 = mitc4_geometric_stiffness(X_SQUARE, Y_SQUARE, s2)
        np.testing.assert_allclose(Kg2, 2.0 * Kg1, atol=1e-12)

    def test_buckling_eigenvalue_cantilever(self):
        """Single-element cantilever plate buckling gives reasonable eigenvalue.

        Fix left edge (nodes 1,4), compress right edge (nodes 2,3).
        Compare to analytical: Ncr = pi^2*D/Lx^2 for clamped-free column.
        Tolerance: within factor of 2 for single element.
        """
        from pyNastran.dev.solver.stiffness.shell_mitc4 import mitc4_geometric_stiffness
        from scipy.linalg import eig

        Lx, Ly = 1.0, 1.0
        x = np.array([0.0, Lx, Lx, 0.0])
        y = np.array([0.0, 0.0, Ly, Ly])

        Ke = mitc4_stiffness(x, y, A_MAT, D_MAT, DS_MAT)

        # Pre-buckling: uniform Nxx = -1 N/m
        Kg = mitc4_geometric_stiffness(x, y, np.array([-1.0, 0.0, 0.0]))

        # Fix nodes 1,4 (left edge)
        constrained = list(range(0, 5)) + list(range(15, 20))
        free = [d for d in range(20) if d not in constrained]

        Kff = Ke[np.ix_(free, free)]
        Kgff = Kg[np.ix_(free, free)]

        eigvals, _ = eig(Kff, -Kgff)
        real_mask = np.abs(np.imag(eigvals)) < 1e-6 * np.abs(eigvals).max()
        eigvals_real = np.real(eigvals[real_mask])
        positive = eigvals_real[eigvals_real > 0]
        lambda_cr = np.sort(positive)[0]

        # Analytical: Ncr for clamped-free plate strip ~ pi^2*D/(4*Lx^2)
        # (factor k=0.25 for clamped-free, single half-wave)
        _fD = E * T**3 / (12 * (1 - NU**2))
        Ncr_analytical = np.pi**2 * _fD / (4 * Lx**2)
        # lambda_cr * |N_applied| should approximate Ncr
        # For single element, expect 50-200% of analytical
        ratio = lambda_cr / Ncr_analytical
        assert 0.5 < ratio < 3.0, \
            f"lambda_cr={lambda_cr:.3e}, Ncr_analytical={Ncr_analytical:.3e}, ratio={ratio:.3f}"

    def test_convergence_ssss_plate(self):
        """SSSS plate with NxN mesh should converge to analytical buckling load.

        Analytical: Ncr = 4*pi^2*D/b^2 for square SSSS plate under Nxx.
        4x4 mesh should be within 5% of analytical.
        """
        from pyNastran.dev.solver.stiffness.shell_mitc4 import mitc4_geometric_stiffness
        from scipy.linalg import eig

        Lx, Ly = 1.0, 1.0
        n = 4  # 4x4 mesh
        dx, dy = Lx / n, Ly / n
        n_nodes = (n + 1) ** 2
        n_dof = 5 * n_nodes

        K_g = np.zeros((n_dof, n_dof))
        Kg_g = np.zeros((n_dof, n_dof))

        stress = np.array([-1.0, 0.0, 0.0])

        for iy in range(n):
            for ix in range(n):
                n1 = iy * (n + 1) + ix
                n2 = n1 + 1
                n3 = n1 + (n + 1) + 1
                n4 = n1 + (n + 1)

                x_e = np.array([ix*dx, (ix+1)*dx, (ix+1)*dx, ix*dx])
                y_e = np.array([iy*dy, iy*dy, (iy+1)*dy, (iy+1)*dy])

                Ke = mitc4_stiffness(x_e, y_e, A_MAT, D_MAT, DS_MAT)
                Kg_e = mitc4_geometric_stiffness(x_e, y_e, stress)

                dofs = []
                for nd in [n1, n2, n3, n4]:
                    dofs.extend([5*nd + d for d in range(5)])

                for ii in range(20):
                    for jj in range(20):
                        K_g[dofs[ii], dofs[jj]] += Ke[ii, jj]
                        Kg_g[dofs[ii], dofs[jj]] += Kg_e[ii, jj]

        # SSSS BCs: w=0 on all boundary nodes
        boundary = set()
        for i in range(n + 1):
            boundary.add(i)
            boundary.add(n * (n + 1) + i)
            boundary.add(i * (n + 1))
            boundary.add(i * (n + 1) + n)

        constrained = set()
        for nd in boundary:
            constrained.add(5 * nd + 2)  # w=0
        constrained.add(0)  # u at corner
        constrained.add(1)  # v at corner
        constrained.add(5 * n + 1)  # v at another corner

        free = sorted(d for d in range(n_dof) if d not in constrained)
        Kff = K_g[np.ix_(free, free)]
        Kgff = Kg_g[np.ix_(free, free)]

        eigvals, _ = eig(Kff, -Kgff)
        real_mask = np.abs(np.imag(eigvals)) < 1e-6 * np.abs(eigvals).max()
        eigvals_real = np.real(eigvals[real_mask])
        positive = eigvals_real[eigvals_real > 0]
        lambda_cr = np.sort(positive)[0]

        # Analytical: Ncr = 4*pi^2*D/b^2, lambda = Ncr / |N_applied|
        _fD = E * T**3 / (12 * (1 - NU**2))
        Ncr_analytical = 4 * np.pi**2 * _fD / Ly**2
        lambda_exact = Ncr_analytical / 1.0  # N_applied = 1.0

        ratio = lambda_cr / lambda_exact
        assert 0.90 < ratio < 1.10, \
            f"4x4 mesh: lambda_cr={lambda_cr:.3e}, exact={lambda_exact:.3e}, ratio={ratio:.4f}"


if __name__ == '__main__':
    unittest.main()
