"""Tests for PG reduction from g-set to a-set.

Verifies:
1. SPC-only: Pa = Pg[free DOFs] (simple partition)
2. MPC: Pn = Pg_n + Gm^T @ Pg_m (force transfer to independent DOFs)
3. OMIT: Pa = Pf_a - Kao @ Koo^{-1} @ Pf_o (static condensation of loads)
4. Full chain: MPC + SPC + OMIT combined
"""
import unittest
import numpy as np
from pyNastran.dev.bdf_vectorized3.loads.reduce_pg import reduce_Pg_to_Pa


class TestReducePg(unittest.TestCase):
    def test_spc_only(self):
        """SPC removal is a simple partition — forces at SPC DOFs are dropped."""
        ndof = 10
        Pg = np.arange(1.0, ndof + 1)

        # Constrain DOFs 0, 1 (first two)
        sset_b = np.zeros(ndof, dtype="bool")
        sset_b[0] = True
        sset_b[1] = True

        Pa = reduce_Pg_to_Pa(Pg, ndof, sset_b=sset_b)

        expected = Pg[~sset_b]  # DOFs 2..9
        assert np.allclose(Pa, expected), f"Pa={Pa}, expected={expected}"
        assert len(Pa) == ndof - 2


    def test_mpc_force_transfer(self):
        """MPC transfers dependent-DOF forces to independent DOFs.

        Setup: 4 DOFs. DOF 0 is dependent on DOFs 1,2 via:
            u_0 = 0.5 * u_1 + 0.5 * u_2
        So Gm = [[0.5, 0.5, 0.0]]  (1 dep DOF, 3 indep DOFs)

        A force of 10.0 at DOF 0 should transfer to:
            Pn_1 += 0.5 * 10 = 5.0
            Pn_2 += 0.5 * 10 = 5.0
        """
        ndof = 4
        Pg = np.array([10.0, 1.0, 2.0, 3.0])

        # DOF 0 is dependent (m-set)
        mset_b = np.zeros(ndof, dtype="bool")
        mset_b[0] = True

        # Gm: (n_m=1, n_n=3) — maps independent DOFs to dependent
        # u_0 = 0.5*u_1 + 0.5*u_2 + 0*u_3
        Gm = np.array([[0.5, 0.5, 0.0]])

        # No SPC, no OMIT
        sset_b = np.zeros(ndof, dtype="bool")

        Pa = reduce_Pg_to_Pa(Pg, ndof, sset_b=sset_b, mset_b=mset_b, Gm=Gm)

        # Pn = Pg_n + Gm^T @ Pg_m
        # Pg_n = [1.0, 2.0, 3.0] (DOFs 1,2,3)
        # Gm^T @ Pg_m = [[0.5],[0.5],[0.0]] * 10 = [5.0, 5.0, 0.0]
        # Pn = [6.0, 7.0, 3.0]
        expected = np.array([6.0, 7.0, 3.0])
        assert np.allclose(Pa, expected), f"Pa={Pa}, expected={expected}"


    def test_omit_load_condensation(self):
        """OMIT condenses interior loads onto boundary via static reduction.

        Setup: 4 free DOFs (no SPC, no MPC). DOFs 2,3 are omitted (interior).
        Stiffness:
            K = [[2, -1, 0, 0],
                 [-1, 2, -1, 0],
                 [0, -1, 2, -1],
                 [0, 0, -1, 1]]

        Kaa = K[0:2, 0:2] = [[2,-1],[-1,2]]
        Koo = K[2:4, 2:4] = [[2,-1],[0, 1]]... wait, let me use proper indexing.

        Actually, let's use a simpler 4-DOF spring chain:
        k=1 between each pair, fixed at left (DOF removed).
        Free DOFs: 0,1,2,3. Omit DOFs 2,3 (interior).

        K_full = [[2,-1,0,0],[-1,2,-1,0],[0,-1,2,-1],[0,0,-1,1]]
        A-set: DOFs 0,1. O-set: DOFs 2,3.
        """
        ndof = 4
        K = np.array([
            [2.0, -1.0, 0.0, 0.0],
            [-1.0, 2.0, -1.0, 0.0],
            [0.0, -1.0, 2.0, -1.0],
            [0.0, 0.0, -1.0, 1.0],
        ])

        # No MPC, no SPC
        sset_b = np.zeros(ndof, dtype="bool")
        mset_b = None

        # OMIT DOFs 2, 3
        oset_b = np.zeros(ndof, dtype="bool")
        oset_b[2] = True
        oset_b[3] = True

        # Extract submatrices
        a_idx = np.array([0, 1])
        o_idx = np.array([2, 3])
        Koo = K[np.ix_(o_idx, o_idx)]
        Koa = K[np.ix_(o_idx, a_idx)]

        # Apply force at interior DOF 3 only
        Pg = np.array([0.0, 0.0, 0.0, 5.0])

        Pa = reduce_Pg_to_Pa(
            Pg, ndof, sset_b=sset_b, mset_b=mset_b, oset_b=oset_b,
            Koo=Koo, Koa=Koa)

        # Manual calculation:
        # Pf_a = [0.0, 0.0], Pf_o = [0.0, 5.0]
        # Koo^{-1} @ Pf_o:
        Koo_inv_Pfo = np.linalg.solve(Koo, np.array([0.0, 5.0]))
        # Kao = Koa^T
        Kao = Koa.T
        expected = np.array([0.0, 0.0]) - Kao @ Koo_inv_Pfo

        assert np.allclose(Pa, expected, atol=1e-12), f"Pa={Pa}, expected={expected}"

        # Verify the result makes physical sense:
        # A force at the tip (DOF 3) should produce a non-zero equivalent
        # force at the boundary DOFs after condensation.
        assert not np.allclose(Pa, 0.0), "Condensed load should be non-zero"


    def test_full_chain_mpc_spc_omit(self):
        """Full reduction chain: MPC + SPC + OMIT.

        6 g-set DOFs:
        - DOF 0: MPC dependent (u_0 = u_1 — rigid link)
        - DOF 5: SPC (constrained)
        - DOFs 3,4: OMIT (interior)
        - DOFs 1,2: a-set (analysis)
        """
        ndof = 6

        # MPC: u_0 = 1.0 * u_1 → Gm = [[1, 0, 0, 0, 0]] (relative to n-set)
        mset_b = np.zeros(ndof, dtype="bool")
        mset_b[0] = True
        # n-set = DOFs 1,2,3,4,5 (5 DOFs)
        # Gm: u_dep = Gm @ u_indep → u_0 = 1*u_1 + 0*u_2 + 0*u_3 + 0*u_4 + 0*u_5
        Gm = np.array([[1.0, 0.0, 0.0, 0.0, 0.0]])

        # SPC: DOF 5 is constrained
        sset_b = np.zeros(ndof, dtype="bool")
        sset_b[5] = True

        # OMIT: DOFs 3,4 are interior
        oset_b = np.zeros(ndof, dtype="bool")
        oset_b[3] = True
        oset_b[4] = True

        # After MPC+SPC: f-set = {1, 2, 3, 4} (4 DOFs)
        # After OMIT: a-set = {1, 2}, o-set = {3, 4}
        # Build a simple stiffness for the f-set partition
        # K_ff for DOFs 1,2,3,4 (spring chain k=1):
        Kff = np.array([
            [2.0, -1.0, 0.0, 0.0],
            [-1.0, 2.0, -1.0, 0.0],
            [0.0, -1.0, 2.0, -1.0],
            [0.0, 0.0, -1.0, 1.0],
        ])
        # a-indices in f-set: [0, 1] (DOFs 1,2)
        # o-indices in f-set: [2, 3] (DOFs 3,4)
        Koo = Kff[2:4, 2:4]
        Koa = Kff[2:4, 0:2]

        # Load: 10.0 at DOF 0 (MPC dependent), 3.0 at DOF 4 (OMIT interior)
        Pg = np.array([10.0, 0.0, 0.0, 0.0, 3.0, 0.0])

        Pa = reduce_Pg_to_Pa(
            Pg, ndof, sset_b=sset_b, mset_b=mset_b, oset_b=oset_b, Gm=Gm,
            Koo=Koo, Koa=Koa)

        # Step by step:
        # 1. MPC: Pn = Pg_n + Gm^T @ Pg_m
        #    Pg_m = [10.0], Pg_n = [0, 0, 0, 3, 0] (DOFs 1,2,3,4,5)
        #    Gm^T @ Pg_m = [[1],[0],[0],[0],[0]] * 10 = [10,0,0,0,0]
        #    Pn = [10, 0, 0, 3, 0]
        Pn_expected = np.array([10.0, 0.0, 0.0, 3.0, 0.0])

        # 2. SPC: remove DOF 5 (last in n-set)
        #    Pf = Pn[:-1] = [10, 0, 0, 3]
        Pf_expected = np.array([10.0, 0.0, 0.0, 3.0])

        # 3. OMIT: Pa = Pf_a - Kao @ Koo^{-1} @ Pf_o
        #    Pf_a = [10, 0], Pf_o = [0, 3]
        Pf_a = Pf_expected[:2]
        Pf_o = Pf_expected[2:]
        Koo_inv_Pfo = np.linalg.solve(Koo, Pf_o)
        Kao = Koa.T
        Pa_expected = Pf_a - Kao @ Koo_inv_Pfo

        assert Pa.shape == (2,), f"Pa shape={Pa.shape}, expected (2,)"
        assert np.allclose(Pa, Pa_expected, atol=1e-12), f"Full chain: Pa={Pa}, expected={Pa_expected}"


    def test_multiple_load_vectors(self):
        """Reduction works for multiple load columns simultaneously."""
        ndof = 6
        sset_b = np.zeros(ndof, dtype="bool")
        sset_b[0] = True  # DOF 0 constrained

        # 3 load cases
        Pg = np.zeros((ndof, 3))
        Pg[1, 0] = 1.0
        Pg[3, 1] = 2.0
        Pg[5, 2] = 3.0

        Pa = reduce_Pg_to_Pa(Pg, ndof, sset_b=sset_b)

        assert Pa.shape == (5, 3), f"Shape={Pa.shape}"
        # Should just be Pg with row 0 removed
        expected = Pg[1:, :]
        assert np.allclose(Pa, expected)


    def test_zero_load_at_spc(self):
        """Forces at SPC DOFs are simply dropped (they become reactions)."""
        ndof = 4
        Pg = np.array([100.0, 0.0, 0.0, 0.0])  # All force at SPC DOF

        sset_b = np.zeros(ndof, dtype="bool")
        sset_b[0] = True

        Pa = reduce_Pg_to_Pa(Pg, ndof, sset_b=sset_b)
        assert np.allclose(Pa, np.zeros(3)), "Force at SPC should be dropped"


if __name__ == "__main__":  # pragma: no cover
    unittest.main()
