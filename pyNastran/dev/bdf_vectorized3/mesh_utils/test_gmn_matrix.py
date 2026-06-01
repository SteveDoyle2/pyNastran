"""Tests for GMN matrix assembly from rigid elements (RBAR, RBAR1, RBE1, RBE2, RBE3, RROD)."""

import numpy as np
from numpy.testing import assert_allclose

from pyNastran.dev.bdf_vectorized3.mesh_utils.gmn_matrix import (
    _expand_dof_string,
    _rigid_body_matrix,
    assemble_gmn,
    assemble_gmn_simple,
    compute_rbar_gmn_row,
    compute_rbe1_gmn_row,
    compute_rbe2_gmn_row,
    compute_rbe3_constraint_matrices,
    compute_rbe3_gmn_row,
    compute_rbe3_thermal_load,
    compute_rrod_gmn_row,
)


# ============================================================================
# Helper utilities
# ============================================================================


def _make_model_with_rbe3(grids, rbe3_cards):
    """Build a minimal vectorized3 BDF model with RBE3s.

    Parameters
    ----------
    grids : list of (nid, x, y, z)
    rbe3_cards : list of (eid, refgrid, refc, wt_cg_groups)
        where wt_cg_groups = [(weight, comp, [nid, ...])]
    """
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

    model = BDF()
    for nid, x, y, z in grids:
        model.add_grid(nid, [x, y, z])
    for eid, refgrid, refc, wt_cg_groups in rbe3_cards:
        weights = []
        comps = []
        gijs = []
        for w, c, nids in wt_cg_groups:
            weights.append(w)
            comps.append(c)
            gijs.append(nids)
        model.add_rbe3(eid, refgrid, refc, weights, comps, gijs)
    model.setup()
    return model


def _make_model_with_rbe2(grids, rbe2_cards):
    """Build a minimal vectorized3 BDF model with RBE2s.

    Parameters
    ----------
    grids : list of (nid, x, y, z)
    rbe2_cards : list of (eid, gn, cm, [gm1, gm2, ...])
    """
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

    model = BDF()
    for nid, x, y, z in grids:
        model.add_grid(nid, [x, y, z])
    for eid, gn, cm, gmi in rbe2_cards:
        model.add_rbe2(eid, gn, cm, gmi)
    model.setup()
    return model


def _make_model_with_rbe2_and_rbe3(grids, rbe2_cards, rbe3_cards):
    """Build a minimal vectorized3 BDF model with both RBE2s and RBE3s."""
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

    model = BDF()
    for nid, x, y, z in grids:
        model.add_grid(nid, [x, y, z])
    for eid, gn, cm, gmi in rbe2_cards:
        model.add_rbe2(eid, gn, cm, gmi)
    for eid, refgrid, refc, wt_cg_groups in rbe3_cards:
        weights = []
        comps = []
        gijs = []
        for w, c, nids in wt_cg_groups:
            weights.append(w)
            comps.append(c)
            gijs.append(nids)
        model.add_rbe3(eid, refgrid, refc, weights, comps, gijs)
    model.setup()
    return model


def _make_model_with_rbe1(grids, rbe1_cards):
    """Build a minimal vectorized3 BDF model with RBE1s.

    Parameters
    ----------
    grids : list of (nid, x, y, z)
    rbe1_cards : list of (eid, Gni, Cni, Gmi, Cmi)
        Gni/Cni = independent (define rigid body), Gmi/Cmi = dependent (UM)
    """
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

    model = BDF()
    for nid, x, y, z in grids:
        model.add_grid(nid, [x, y, z])
    for eid, Gni, Cni, Gmi, Cmi in rbe1_cards:
        model.add_rbe1(eid, Gni, Cni, Gmi, Cmi)
    model.setup()
    return model


def _make_model_with_rbar(grids, rbar_cards):
    """Build a minimal vectorized3 BDF model with RBARs.

    Parameters
    ----------
    grids : list of (nid, x, y, z)
    rbar_cards : list of (eid, [ga, gb], cna, cnb, cma, cmb)
    """
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

    model = BDF()
    for nid, x, y, z in grids:
        model.add_grid(nid, [x, y, z])
    for eid, nids, cna, cnb, cma, cmb in rbar_cards:
        model.add_rbar(eid, nids, cna, cnb, cma, cmb)
    model.setup()
    return model


def _make_model_with_rbar1(grids, rbar1_cards):
    """Build a minimal vectorized3 BDF model with RBAR1s.

    Parameters
    ----------
    grids : list of (nid, x, y, z)
    rbar1_cards : list of (eid, [ga, gb], cb)
    """
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

    model = BDF()
    for nid, x, y, z in grids:
        model.add_grid(nid, [x, y, z])
    for eid, nids, cb in rbar1_cards:
        model.add_rbar1(eid, nids, cb)
    model.setup()
    return model


def _make_model_with_rrod(grids, rrod_cards):
    """Build a minimal vectorized3 BDF model with RRODs.

    Parameters
    ----------
    grids : list of (nid, x, y, z)
    rrod_cards : list of (eid, [ga, gb], cma, cmb)
    """
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

    model = BDF()
    for nid, x, y, z in grids:
        model.add_grid(nid, [x, y, z])
    for eid, nids, cma, cmb in rrod_cards:
        model.add_rrod(eid, nids, cma, cmb)
    model.setup()
    return model


def _make_model_all_types(
    grids,
    rbe1_cards=None,
    rbe2_cards=None,
    rbe3_cards=None,
    rbar_cards=None,
    rbar1_cards=None,
    rrod_cards=None,
):
    """Build a model with any combination of rigid element types."""
    from pyNastran.dev.bdf_vectorized3.bdf import BDF

    model = BDF()
    for nid, x, y, z in grids:
        model.add_grid(nid, [x, y, z])
    if rbe1_cards:
        for eid, Gni, Cni, Gmi, Cmi in rbe1_cards:
            model.add_rbe1(eid, Gni, Cni, Gmi, Cmi)
    if rbe2_cards:
        for eid, gn, cm, gmi in rbe2_cards:
            model.add_rbe2(eid, gn, cm, gmi)
    if rbe3_cards:
        for eid, refgrid, refc, wt_cg_groups in rbe3_cards:
            weights = []
            comps = []
            gijs = []
            for w, c, nids in wt_cg_groups:
                weights.append(w)
                comps.append(c)
                gijs.append(nids)
            model.add_rbe3(eid, refgrid, refc, weights, comps, gijs)
    if rbar_cards:
        for eid, nids, cna, cnb, cma, cmb in rbar_cards:
            model.add_rbar(eid, nids, cna, cnb, cma, cmb)
    if rbar1_cards:
        for eid, nids, cb in rbar1_cards:
            model.add_rbar1(eid, nids, cb)
    if rrod_cards:
        for eid, nids, cma, cmb in rrod_cards:
            model.add_rrod(eid, nids, cma, cmb)
    model.setup()
    return model


# ============================================================================
# Tests: utility functions
# ============================================================================


def test_expand_dof_string():
    """_expand_dof_string converts component integers to 0-based lists."""
    assert _expand_dof_string(123) == [0, 1, 2]
    assert _expand_dof_string(456) == [3, 4, 5]
    assert _expand_dof_string(123456) == [0, 1, 2, 3, 4, 5]
    assert _expand_dof_string(1) == [0]
    assert _expand_dof_string(36) == [2, 5]
    assert _expand_dof_string("123") == [0, 1, 2]


def test_rigid_body_matrix_zero_offset():
    """With zero offset, RB matrix is identity."""
    T = _rigid_body_matrix(np.zeros(3))
    assert_allclose(T, np.eye(6), atol=1e-15)


def test_rigid_body_matrix_offset():
    """RB matrix with non-zero offset has coupling terms."""
    dx = np.array([1.0, 2.0, 3.0])
    T = _rigid_body_matrix(dx)
    # Check diagonal is 1
    assert_allclose(np.diag(T), np.ones(6), atol=1e-15)
    # Translation-rotation coupling: T[0,4]=z, T[0,5]=-y
    assert_allclose(T[0, 4], 3.0, atol=1e-15)
    assert_allclose(T[0, 5], -2.0, atol=1e-15)
    # T[1,3]=-z, T[1,5]=x
    assert_allclose(T[1, 3], -3.0, atol=1e-15)
    assert_allclose(T[1, 5], 1.0, atol=1e-15)
    # T[2,3]=y, T[2,4]=-x
    assert_allclose(T[2, 3], 2.0, atol=1e-15)
    assert_allclose(T[2, 4], -1.0, atol=1e-15)


# ============================================================================
# Tests: RBE3 GMN
# ============================================================================


def test_rbe3_gmn_collinear_2node_translations():
    """Two collinear independent nodes, equal weights, translations only.

    Dependent node at midpoint => GMN should average the translations.
    """
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array(
        [
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
    )
    dep_dofs = [0, 1, 2]  # T1, T2, T3
    ind_dofs_per_node = [[0, 1, 2], [0, 1, 2]]
    weights = np.array([1.0, 1.0])

    gmn, _ = compute_rbe3_gmn_row(xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights)

    # Shape: (3 dep dofs, 6 ind dofs total)
    assert gmn.shape == (3, 6)
    # Each translation should be 0.5 from each node
    expected = np.array(
        [
            [0.5, 0, 0, 0.5, 0, 0],
            [0, 0.5, 0, 0, 0.5, 0],
            [0, 0, 0.5, 0, 0, 0.5],
        ]
    )
    assert_allclose(gmn, expected, atol=1e-12)


def test_rbe3_gmn_single_node():
    """Single independent node: GMN = identity (selected DOFs)."""
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array([[0.0, 0.0, 0.0]])
    dep_dofs = [0, 1, 2, 3, 4, 5]
    ind_dofs_per_node = [[0, 1, 2, 3, 4, 5]]
    weights = np.array([1.0])

    gmn, _ = compute_rbe3_gmn_row(xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights)

    assert gmn.shape == (6, 6)
    assert_allclose(gmn, np.eye(6), atol=1e-12)


def test_rbe3_gmn_single_node_offset():
    """Single independent node with offset: GMN = rigid-body transform."""
    xyz_dep = np.array([0.0, 0.0, 0.0])
    xyz_ind = np.array([1.0, 0.0, 0.0])
    ind_nodes_xyz = np.array([xyz_ind])
    dep_dofs = [0, 1, 2, 3, 4, 5]
    ind_dofs_per_node = [[0, 1, 2, 3, 4, 5]]
    weights = np.array([1.0])

    gmn, _ = compute_rbe3_gmn_row(xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights)

    # The GMN should be the inverse rigid-body relation
    # dep motion = GMN * ind motion
    # For single node with full DOFs, GMN = T^{-1} where T maps dep->ind
    # Actually for a single node, it inverts to the RB matrix from ind to dep
    dx = xyz_ind - xyz_dep  # (1, 0, 0)
    T = _rigid_body_matrix(dx)  # maps dep 6DOF -> ind 6DOF
    # GMN = inv(T^T W T) @ T^T W = inv(T) since W=I and T is invertible
    expected = np.linalg.inv(T)
    assert_allclose(gmn, expected, atol=1e-12)


def test_rbe3_gmn_unequal_weights():
    """Two nodes with unequal weights: weighted interpolation.

    For collinear nodes along x with translation-only DOFs, the T1
    weighting follows w1/(w1+w2) directly. T2/T3 have rotation coupling
    through the rank-deficient normal equations.
    """
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array(
        [
            [-1.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
    )
    dep_dofs = [0, 1, 2]
    ind_dofs_per_node = [[0, 1, 2], [0, 1, 2]]
    weights = np.array([1.0, 3.0])

    gmn, _ = compute_rbe3_gmn_row(xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights)

    assert gmn.shape == (3, 6)
    # T1 (x-dir): no rotation coupling for collinear x-offset
    assert_allclose(gmn[0, 0], 0.25, atol=1e-12)
    assert_allclose(gmn[0, 3], 0.75, atol=1e-12)
    # Row sums must equal 1 for each dep DOF (rigid body translation preserved)
    for i in range(3):
        assert_allclose(gmn[i, i] + gmn[i, i + 3], 1.0, atol=1e-12)


def test_rbe3_gmn_symmetric_4nodes():
    """Four symmetric nodes around dependent => equal weights."""
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ]
    )
    dep_dofs = [0, 1, 2]
    ind_dofs_per_node = [[0, 1, 2]] * 4
    weights = np.ones(4)

    gmn, _ = compute_rbe3_gmn_row(xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights)

    # By symmetry, each translation DOF gets 0.25 from each node
    assert gmn.shape == (3, 12)
    for i_dep in range(3):
        for i_node in range(4):
            col = i_node * 3 + i_dep
            assert_allclose(gmn[i_dep, col], 0.25, atol=1e-12)


def test_rbe3_constraint_matrices():
    """compute_rbe3_constraint_matrices returns RMM=I, RMN=-GMN."""
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
        ]
    )
    dep_dofs = [0, 1, 2, 3, 4, 5]
    ind_dofs_per_node = [[0, 1, 2, 3, 4, 5]] * 2
    weights = np.ones(2)

    RMM, RMN, GMN = compute_rbe3_constraint_matrices(
        xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights
    )

    assert_allclose(RMM, np.eye(6), atol=1e-15)
    assert_allclose(RMN, -GMN, atol=1e-15)


def test_rbe3_assemble_gmn_simple():
    """assemble_gmn_simple interface produces correct shape and values."""
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
        ]
    )

    gmn = assemble_gmn_simple(xyz_dep, 123456, ind_nodes_xyz, 123456)
    assert gmn.shape == (6, 12)

    # Translations should average (equal weights, symmetric about dep)
    assert_allclose(gmn[0, 0], 0.5, atol=1e-12)
    assert_allclose(gmn[0, 6], 0.5, atol=1e-12)


def test_rbe3_thermal_load_zero_alpha():
    """Zero alpha or zero dT => zero thermal displacement."""
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array([[1.0, 0.0, 0.0]])
    dep_dofs = [0, 1, 2]
    ind_dofs_per_node = [[0, 1, 2]]
    weights = np.array([1.0])

    gmn, _ = compute_rbe3_gmn_row(xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights)

    u_th = compute_rbe3_thermal_load(
        gmn,
        xyz_dep,
        dep_dofs,
        ind_nodes_xyz,
        ind_dofs_per_node,
        alpha=0.0,
        tref=70.0,
        temperature=200.0,
    )
    assert_allclose(u_th, np.zeros(3), atol=1e-15)

    u_th2 = compute_rbe3_thermal_load(
        gmn,
        xyz_dep,
        dep_dofs,
        ind_nodes_xyz,
        ind_dofs_per_node,
        alpha=1e-5,
        tref=70.0,
        temperature=70.0,
    )
    assert_allclose(u_th2, np.zeros(3), atol=1e-15)


def test_rbe3_thermal_load_nonzero():
    """Non-zero alpha and dT produce expected thermal displacement."""
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array([[2.0, 0.0, 0.0]])
    dep_dofs = [0, 1, 2]
    ind_dofs_per_node = [[0, 1, 2]]
    weights = np.array([1.0])

    gmn, _ = compute_rbe3_gmn_row(xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights)

    alpha = 1e-5
    tref = 70.0
    temperature = 170.0
    dT = temperature - tref  # 100

    u_th = compute_rbe3_thermal_load(
        gmn,
        xyz_dep,
        dep_dofs,
        ind_nodes_xyz,
        ind_dofs_per_node,
        alpha=alpha,
        tref=tref,
        temperature=temperature,
    )
    # u_thermal[T1] = alpha * dT * dx[0] = 1e-5 * 100 * 2.0 = 0.002
    # For single node with comp 123 at dep, gmn = I (3x3),
    # so u_thermal_dep = gmn @ u_thermal = u_thermal
    assert_allclose(u_th[0], alpha * dT * 2.0, atol=1e-15)
    assert_allclose(u_th[1], 0.0, atol=1e-15)
    assert_allclose(u_th[2], 0.0, atol=1e-15)


def test_rbe3_assemble_gmn_simple_with_thermal():
    """assemble_gmn_simple returns (gmn, u_thermal) tuple when T given."""
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array([[1.0, 0.0, 0.0]])

    result = assemble_gmn_simple(
        xyz_dep, 123, ind_nodes_xyz, 123, alpha=1e-5, tref=70.0, temperature=170.0
    )
    assert isinstance(result, tuple)
    gmn, u_th = result
    assert gmn.shape == (3, 3)
    assert u_th.shape == (3,)
    assert_allclose(u_th[0], 1e-5 * 100.0 * 1.0, atol=1e-15)


def test_rbe3_rigid_body_preservation():
    """RBE3 GMN preserves rigid-body motion: if all independent nodes move
    rigidly, the dependent node moves the same way.
    """
    xyz_dep = np.array([0.0, 0.0, 0.0])
    ind_nodes_xyz = np.array(
        [
            [1.0, 0.0, 0.0],
            [-1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, -1.0, 0.0],
        ]
    )
    dep_dofs = [0, 1, 2, 3, 4, 5]
    ind_dofs_per_node = [[0, 1, 2, 3, 4, 5]] * 4
    weights = np.ones(4)

    gmn, _ = compute_rbe3_gmn_row(xyz_dep, dep_dofs, ind_nodes_xyz, ind_dofs_per_node, weights)

    # Apply a pure translation + rotation
    # Translation: (1, 2, 3), Rotation about origin: (0.01, 0.02, 0.03)
    u_dep_ref = np.array([1.0, 2.0, 3.0, 0.01, 0.02, 0.03])
    rx, ry, rz = 0.01, 0.02, 0.03
    tx, ty, tz = 1.0, 2.0, 3.0

    u_ind = np.zeros(24)
    for i_node in range(4):
        dx = ind_nodes_xyz[i_node]
        # Rigid body: u = T + omega x r
        u_ind[i_node * 6 + 0] = tx + ry * dx[2] - rz * dx[1]
        u_ind[i_node * 6 + 1] = ty - rx * dx[2] + rz * dx[0]
        u_ind[i_node * 6 + 2] = tz + rx * dx[1] - ry * dx[0]
        u_ind[i_node * 6 + 3] = rx
        u_ind[i_node * 6 + 4] = ry
        u_ind[i_node * 6 + 5] = rz

    u_dep_computed = gmn @ u_ind
    assert_allclose(u_dep_computed, u_dep_ref, atol=1e-12)


def test_rbe3_global_assembly():
    """Test assemble_gmn with a model containing one RBE3."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
        (3, -1.0, 0.0, 0.0),
    ]
    rbe3_cards = [
        (100, 1, "123456", [(1.0, "123", [2, 3])]),
    ]
    model = _make_model_with_rbe3(grids, rbe3_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    # m-set: 6 DOFs of node 1
    assert len(m_set) == 6
    for dof in range(1, 7):
        assert (1, dof) in m_set

    # GMN shape: (6 m-set DOFs, 18 total DOFs)
    assert GMN.shape == (6, 18)

    # Verify translation average for T1
    m_idx = m_set[(1, 1)]
    col_2_t1 = dof_map[(2, 1)]
    col_3_t1 = dof_map[(3, 1)]
    assert_allclose(GMN[m_idx, col_2_t1], 0.5, atol=1e-12)
    assert_allclose(GMN[m_idx, col_3_t1], 0.5, atol=1e-12)


# ============================================================================
# Tests: RBE2 GMN
# ============================================================================


def test_rbe2_gmn_coincident():
    """RBE2 with coincident nodes: GMN is identity for constrained DOFs."""
    xyz_ind = np.array([0.0, 0.0, 0.0])
    xyz_dep = np.array([0.0, 0.0, 0.0])
    dep_dofs = [0, 1, 2, 3, 4, 5]

    gmn = compute_rbe2_gmn_row(xyz_ind, xyz_dep, dep_dofs)

    assert gmn.shape == (6, 6)
    assert_allclose(gmn, np.eye(6), atol=1e-15)


def test_rbe2_gmn_offset_x():
    """RBE2 with x-offset: dependent translations couple to rotations."""
    xyz_ind = np.array([0.0, 0.0, 0.0])
    xyz_dep = np.array([1.0, 0.0, 0.0])
    dep_dofs = [0, 1, 2, 3, 4, 5]

    gmn = compute_rbe2_gmn_row(xyz_ind, xyz_dep, dep_dofs)

    # offset = xyz_dep - xyz_ind = (1, 0, 0)
    # T = rigid_body_matrix([1, 0, 0])
    expected = _rigid_body_matrix(np.array([1.0, 0.0, 0.0]))
    assert_allclose(gmn, expected, atol=1e-15)

    # Check specific coupling: T[1,5] = x = 1.0 (T2 from R3)
    assert_allclose(gmn[1, 5], 1.0, atol=1e-15)
    # T[2,4] = -x = -1.0 (T3 from R2)
    assert_allclose(gmn[2, 4], -1.0, atol=1e-15)


def test_rbe2_gmn_offset_xyz():
    """RBE2 with general offset."""
    xyz_ind = np.array([1.0, 2.0, 3.0])
    xyz_dep = np.array([4.0, 5.0, 6.0])
    dep_dofs = [0, 1, 2, 3, 4, 5]

    gmn = compute_rbe2_gmn_row(xyz_ind, xyz_dep, dep_dofs)

    offset = xyz_dep - xyz_ind  # (3, 3, 3)
    expected = _rigid_body_matrix(offset)
    assert_allclose(gmn, expected, atol=1e-15)


def test_rbe2_gmn_partial_dofs():
    """RBE2 with only translational DOFs constrained (CM=123)."""
    xyz_ind = np.array([0.0, 0.0, 0.0])
    xyz_dep = np.array([1.0, 0.0, 0.0])
    dep_dofs = [0, 1, 2]

    gmn = compute_rbe2_gmn_row(xyz_ind, xyz_dep, dep_dofs)

    assert gmn.shape == (3, 6)
    expected_full = _rigid_body_matrix(np.array([1.0, 0.0, 0.0]))
    assert_allclose(gmn, expected_full[:3, :], atol=1e-15)


def test_rbe2_gmn_rotation_only():
    """RBE2 with only rotational DOFs constrained (CM=456)."""
    xyz_ind = np.array([0.0, 0.0, 0.0])
    xyz_dep = np.array([1.0, 2.0, 3.0])
    dep_dofs = [3, 4, 5]

    gmn = compute_rbe2_gmn_row(xyz_ind, xyz_dep, dep_dofs)

    assert gmn.shape == (3, 6)
    # Rotational rows have no coupling to translations in the RB matrix
    # Rows 3,4,5 of the RB matrix are just identity for the rotation block
    expected = np.array(
        [
            [0, 0, 0, 1, 0, 0],
            [0, 0, 0, 0, 1, 0],
            [0, 0, 0, 0, 0, 1],
        ],
        dtype=float,
    )
    assert_allclose(gmn, expected, atol=1e-15)


def test_rbe2_rigid_body_preservation():
    """RBE2: rigid body motion at independent node reproduces exactly
    at the dependent node.
    """
    xyz_ind = np.array([0.0, 0.0, 0.0])
    xyz_dep = np.array([2.0, 1.0, 0.5])
    dep_dofs = [0, 1, 2, 3, 4, 5]

    gmn = compute_rbe2_gmn_row(xyz_ind, xyz_dep, dep_dofs)

    # Apply rigid body motion at independent node
    tx, ty, tz = 0.1, -0.2, 0.3
    rx, ry, rz = 0.01, -0.02, 0.03
    u_ind = np.array([tx, ty, tz, rx, ry, rz])

    u_dep = gmn @ u_ind

    # Expected dependent motion: translation + omega x offset
    offset = xyz_dep - xyz_ind  # (2, 1, 0.5)
    expected_tx = tx + ry * offset[2] - rz * offset[1]
    expected_ty = ty - rx * offset[2] + rz * offset[0]
    expected_tz = tz + rx * offset[1] - ry * offset[0]
    expected = np.array([expected_tx, expected_ty, expected_tz, rx, ry, rz])

    assert_allclose(u_dep, expected, atol=1e-14)


def test_rbe2_multiple_dependent_nodes():
    """RBE2 with multiple dependent nodes: each gets its own GMN rows."""
    xyz_ind = np.array([0.0, 0.0, 0.0])
    dep_nodes_xyz = np.array(
        [
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ]
    )
    dep_dofs = [0, 1, 2, 3, 4, 5]

    # Apply same rigid body motion to all dependent nodes
    u_ind = np.array([0.5, -0.3, 0.1, 0.01, -0.02, 0.005])

    for xyz_dep in dep_nodes_xyz:
        gmn = compute_rbe2_gmn_row(xyz_ind, xyz_dep, dep_dofs)
        u_dep = gmn @ u_ind

        offset = xyz_dep - xyz_ind
        rx, ry, rz = u_ind[3], u_ind[4], u_ind[5]
        expected_t = u_ind[:3] + np.cross(np.array([rx, ry, rz]), offset)
        expected = np.concatenate([expected_t, u_ind[3:]])
        assert_allclose(u_dep, expected, atol=1e-14)


def test_rbe2_global_assembly_single():
    """Test assemble_gmn with a single RBE2 element."""
    grids = [
        (1, 0.0, 0.0, 0.0),  # independent
        (2, 1.0, 0.0, 0.0),  # dependent
        (3, 0.0, 1.0, 0.0),  # dependent
    ]
    rbe2_cards = [
        (100, 1, "123456", [2, 3]),
    ]
    model = _make_model_with_rbe2(grids, rbe2_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    # 12 dependent DOFs (2 nodes x 6 DOFs)
    assert len(m_set) == 12
    assert GMN.shape == (12, 18)

    # Check node 2 T1 DOF: maps from node 1's DOFs
    m_idx = m_set[(2, 1)]
    col_1_t1 = dof_map[(1, 1)]
    assert_allclose(GMN[m_idx, col_1_t1], 1.0, atol=1e-15)

    # Check node 2 T2: has coupling from node 1 R3 (offset x=1)
    m_t2 = m_set[(2, 2)]
    col_1_r3 = dof_map[(1, 6)]
    assert_allclose(GMN[m_t2, col_1_r3], 1.0, atol=1e-15)

    # Check node 3 T1: has coupling from node 1 R3 (offset y=1)
    m_3_t1 = m_set[(3, 1)]
    col_1_r3 = dof_map[(1, 6)]
    # T[0,5] = -y = -1.0
    assert_allclose(GMN[m_3_t1, col_1_r3], -1.0, atol=1e-15)


def test_rbe2_global_assembly_rigid_motion():
    """Verify rigid-body preservation through the global GMN assembly."""
    grids = [
        (10, 0.0, 0.0, 0.0),
        (20, 3.0, 0.0, 0.0),
        (30, 0.0, 4.0, 0.0),
    ]
    rbe2_cards = [
        (1, 10, "123456", [20, 30]),
    ]
    model = _make_model_with_rbe2(grids, rbe2_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    # Build full displacement vector with rigid body motion at node 10
    ndof = 18
    u_full = np.zeros(ndof)
    tx, ty, tz = 0.5, -0.3, 0.1
    rx, ry, rz = 0.01, -0.005, 0.02

    # Set independent node DOFs
    for i, val in enumerate([tx, ty, tz, rx, ry, rz]):
        u_full[dof_map[(10, i + 1)]] = val

    # Compute dependent displacements
    u_m = GMN @ u_full

    # Check node 20 (offset = [3, 0, 0])
    offset_20 = np.array([3.0, 0.0, 0.0])
    exp_t20 = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset_20)
    for i_dof in range(3):
        m_idx = m_set[(20, i_dof + 1)]
        assert_allclose(u_m[m_idx], exp_t20[i_dof], atol=1e-14)
    for i_dof in range(3):
        m_idx = m_set[(20, i_dof + 4)]
        assert_allclose(u_m[m_idx], [rx, ry, rz][i_dof], atol=1e-14)

    # Check node 30 (offset = [0, 4, 0])
    offset_30 = np.array([0.0, 4.0, 0.0])
    exp_t30 = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset_30)
    for i_dof in range(3):
        m_idx = m_set[(30, i_dof + 1)]
        assert_allclose(u_m[m_idx], exp_t30[i_dof], atol=1e-14)


# ============================================================================
# Tests: RBE1 GMN
# ============================================================================


def test_rbe1_gmn_single_node_6dof():
    """RBE1 with one GN node providing all 6 DOFs: GMN = rigid-body transform."""
    ind_nodes_xyz = np.array([[0.0, 0.0, 0.0]])
    ind_dofs_per_node = [[0, 1, 2, 3, 4, 5]]
    dep_nodes_xyz = np.array([[1.0, 0.0, 0.0]])
    dep_dofs_per_node = [[0, 1, 2, 3, 4, 5]]

    gmn = compute_rbe1_gmn_row(ind_nodes_xyz, ind_dofs_per_node, dep_nodes_xyz, dep_dofs_per_node)

    # Same as RBE2: rigid transform from ind to dep
    expected = _rigid_body_matrix(np.array([1.0, 0.0, 0.0]))
    assert_allclose(gmn, expected, atol=1e-12)


def test_rbe1_gmn_two_nodes_split_dofs():
    """RBE1 with two GN nodes each providing 3 DOFs.

    Node 1 at origin provides T1,T2,T3.
    Node 2 at (1,0,0) provides R1,R2,R3.
    Dependent node at (0,1,0) gets all 6 DOFs.
    """
    ind_nodes_xyz = np.array(
        [
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
        ]
    )
    ind_dofs_per_node = [[0, 1, 2], [3, 4, 5]]
    dep_nodes_xyz = np.array([[0.0, 1.0, 0.0]])
    dep_dofs_per_node = [[0, 1, 2, 3, 4, 5]]

    gmn = compute_rbe1_gmn_row(ind_nodes_xyz, ind_dofs_per_node, dep_nodes_xyz, dep_dofs_per_node)

    # 6 dep DOFs, 6 ind DOFs (3+3)
    assert gmn.shape == (6, 6)

    # Verify with rigid body motion
    tx, ty, tz = 0.1, -0.2, 0.3
    rx, ry, rz = 0.01, -0.02, 0.03

    # Independent DOFs: T1,T2,T3 at node1, R1,R2,R3 at node2
    u_ind = np.array([tx, ty, tz, rx, ry, rz])

    u_dep = gmn @ u_ind

    # Expected: dependent node at (0,1,0)
    # offset from ref (centroid of ind nodes = (0.5,0,0)) to dep = (-0.5, 1, 0)
    # But the actual motion is rigid body:
    # u_dep_T = u_T + omega x (xyz_dep - xyz_ref_rigid_body)
    # Since the rigid body is defined by TXY at origin + R at (1,0,0),
    # the rotation is the same everywhere, so:
    dep_offset_from_origin = np.array([0.0, 1.0, 0.0])
    expected_t = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], dep_offset_from_origin)
    expected_r = np.array([rx, ry, rz])
    expected = np.concatenate([expected_t, expected_r])
    assert_allclose(u_dep, expected, atol=1e-12)


def test_rbe1_gmn_rigid_body_preservation():
    """RBE1: rigid body motion at independent DOFs reproduces at dependent."""
    ind_nodes_xyz = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
        ]
    )
    ind_dofs_per_node = [[0, 1, 2], [0, 1, 2]]
    dep_nodes_xyz = np.array(
        [
            [1.0, 1.0, 0.0],
            [1.0, -1.0, 0.0],
        ]
    )
    dep_dofs_per_node = [[0, 1, 2], [0, 1, 2]]

    gmn = compute_rbe1_gmn_row(ind_nodes_xyz, ind_dofs_per_node, dep_nodes_xyz, dep_dofs_per_node)

    # Apply rigid body: translation + small rotation about z
    tx, ty, tz = 1.0, 2.0, 3.0
    rx, ry, rz = 0.0, 0.0, 0.01

    # Independent DOFs: T1,T2,T3 at node1 and T1,T2,T3 at node2
    u_ind = np.zeros(6)
    for i_node, xyz in enumerate(ind_nodes_xyz):
        offset = xyz
        u_t = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset)
        u_ind[i_node * 3 : (i_node + 1) * 3] = u_t

    u_dep = gmn @ u_ind

    # Check each dependent node
    for i_node, xyz in enumerate(dep_nodes_xyz):
        offset = xyz
        expected_t = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset)
        assert_allclose(u_dep[i_node * 3 : (i_node + 1) * 3], expected_t, atol=1e-12)


def test_rbe1_gmn_partial_dependent():
    """RBE1 with dependent node having only translational DOFs."""
    ind_nodes_xyz = np.array([[0.0, 0.0, 0.0]])
    ind_dofs_per_node = [[0, 1, 2, 3, 4, 5]]
    dep_nodes_xyz = np.array([[1.0, 0.0, 0.0]])
    dep_dofs_per_node = [[0, 1, 2]]  # only translations dependent

    gmn = compute_rbe1_gmn_row(ind_nodes_xyz, ind_dofs_per_node, dep_nodes_xyz, dep_dofs_per_node)

    assert gmn.shape == (3, 6)

    # Verify rigid body
    u_ind = np.array([0.1, -0.2, 0.3, 0.01, -0.02, 0.03])
    u_dep = gmn @ u_ind
    offset = np.array([1.0, 0.0, 0.0])
    expected_t = u_ind[:3] + np.cross(u_ind[3:], offset)
    assert_allclose(u_dep, expected_t, atol=1e-12)


def test_rbe1_global_assembly():
    """Test assemble_gmn with a model containing one RBE1."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 2.0, 0.0, 0.0),
        (3, 1.0, 1.0, 0.0),
    ]
    # GN nodes: 1 (123456) defines the rigid body
    # GM node: 3 (123456) is dependent
    rbe1_cards = [
        (100, [1, 2], ["123", "456"], [3], ["123456"]),
    ]
    model = _make_model_with_rbe1(grids, rbe1_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    # 6 dependent DOFs from node 3
    assert len(m_set) == 6
    for dof in range(1, 7):
        assert (3, dof) in m_set

    # Total DOFs: 3 nodes x 6 = 18
    assert GMN.shape == (6, 18)

    # Verify rigid body through global assembly
    ndof = 18
    u_full = np.zeros(ndof)
    tx, ty, tz = 0.5, -0.3, 0.1
    rx, ry, rz = 0.01, -0.005, 0.02

    # Set node 1 DOFs (T1,T2,T3)
    for i, val in enumerate([tx, ty, tz]):
        u_full[dof_map[(1, i + 1)]] = val
    # Set node 2 DOFs (R1,R2,R3) — but also need T at node 2 from rigid body
    offset_2 = np.array([2.0, 0.0, 0.0])
    u_t2 = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset_2)
    for i, val in enumerate(u_t2):
        u_full[dof_map[(2, i + 1)]] = val
    for i, val in enumerate([rx, ry, rz]):
        u_full[dof_map[(2, i + 4)]] = val

    u_m = GMN @ u_full

    # Node 3 at (1,1,0): expected rigid body motion
    offset_3 = np.array([1.0, 1.0, 0.0])
    exp_t3 = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset_3)
    for i_dof in range(3):
        m_idx = m_set[(3, i_dof + 1)]
        assert_allclose(u_m[m_idx], exp_t3[i_dof], atol=1e-12)
    for i_dof in range(3):
        m_idx = m_set[(3, i_dof + 4)]
        assert_allclose(u_m[m_idx], [rx, ry, rz][i_dof], atol=1e-12)


# ============================================================================
# Tests: RBAR GMN
# ============================================================================


def test_rbar_gmn_all_ind_at_a():
    """RBAR with all 6 independent DOFs at node A (like RBE2)."""
    xyz_a = np.array([0.0, 0.0, 0.0])
    xyz_b = np.array([1.0, 0.0, 0.0])
    ind_dofs_a = [0, 1, 2, 3, 4, 5]
    ind_dofs_b = []
    dep_dofs_a = []
    dep_dofs_b = [0, 1, 2, 3, 4, 5]

    gmn = compute_rbar_gmn_row(xyz_a, xyz_b, ind_dofs_a, ind_dofs_b, dep_dofs_a, dep_dofs_b)

    # Same as RBE2: dependent at B, independent at A
    expected = _rigid_body_matrix(np.array([1.0, 0.0, 0.0]))
    assert_allclose(gmn, expected, atol=1e-12)


def test_rbar_gmn_split_dofs():
    """RBAR with DOFs split between nodes A and B.

    CNA=123 (translations at A), CNB=456 (rotations at B).
    CMA=456 (rotations at A), CMB=123 (translations at B).
    """
    xyz_a = np.array([0.0, 0.0, 0.0])
    xyz_b = np.array([2.0, 0.0, 0.0])
    ind_dofs_a = [0, 1, 2]
    ind_dofs_b = [3, 4, 5]
    dep_dofs_a = [3, 4, 5]
    dep_dofs_b = [0, 1, 2]

    gmn = compute_rbar_gmn_row(xyz_a, xyz_b, ind_dofs_a, ind_dofs_b, dep_dofs_a, dep_dofs_b)

    # 6 dep DOFs (3 at A + 3 at B), 6 ind DOFs (3 at A + 3 at B)
    assert gmn.shape == (6, 6)

    # Verify with rigid body motion
    tx, ty, tz = 0.5, -0.3, 0.1
    rx, ry, rz = 0.01, -0.005, 0.02

    # Independent: T at A, R at B (rotations are same everywhere for rigid body)
    u_ind = np.array([tx, ty, tz, rx, ry, rz])

    u_dep = gmn @ u_ind

    # Dependent: R at A (should be rx,ry,rz), T at B
    offset_b = np.array([2.0, 0.0, 0.0])
    exp_t_b = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset_b)

    # dep_dofs_a = R at A (rows 0-2), dep_dofs_b = T at B (rows 3-5)
    assert_allclose(u_dep[:3], [rx, ry, rz], atol=1e-12)
    assert_allclose(u_dep[3:], exp_t_b, atol=1e-12)


def test_rbar_gmn_rigid_body_preservation():
    """RBAR: rigid body motion is preserved."""
    xyz_a = np.array([0.0, 0.0, 0.0])
    xyz_b = np.array([3.0, 1.0, -1.0])
    ind_dofs_a = [0, 1, 2, 3, 4, 5]
    ind_dofs_b = []
    dep_dofs_a = []
    dep_dofs_b = [0, 1, 2, 3, 4, 5]

    gmn = compute_rbar_gmn_row(xyz_a, xyz_b, ind_dofs_a, ind_dofs_b, dep_dofs_a, dep_dofs_b)

    u_ind = np.array([0.1, -0.2, 0.3, 0.01, -0.02, 0.03])
    u_dep = gmn @ u_ind

    offset = xyz_b - xyz_a
    expected_t = u_ind[:3] + np.cross(u_ind[3:], offset)
    expected = np.concatenate([expected_t, u_ind[3:]])
    assert_allclose(u_dep, expected, atol=1e-12)


def test_rbar_gmn_all_ind_at_b():
    """RBAR with all 6 independent DOFs at node B."""
    xyz_a = np.array([0.0, 0.0, 0.0])
    xyz_b = np.array([0.0, 2.0, 0.0])
    ind_dofs_a = []
    ind_dofs_b = [0, 1, 2, 3, 4, 5]
    dep_dofs_a = [0, 1, 2, 3, 4, 5]
    dep_dofs_b = []

    gmn = compute_rbar_gmn_row(xyz_a, xyz_b, ind_dofs_a, ind_dofs_b, dep_dofs_a, dep_dofs_b)

    # offset from B to A = A - B = (0, -2, 0)
    u_ind = np.array([0.1, -0.2, 0.3, 0.01, -0.02, 0.03])
    u_dep = gmn @ u_ind

    offset_a_from_b = xyz_a - xyz_b  # (0, -2, 0)
    expected_t = u_ind[:3] + np.cross(u_ind[3:], offset_a_from_b)
    expected = np.concatenate([expected_t, u_ind[3:]])
    assert_allclose(u_dep, expected, atol=1e-12)


def test_rbar_global_assembly():
    """Test assemble_gmn with an RBAR element."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
    ]
    # All independent at node 1, all dependent at node 2
    rbar_cards = [
        (100, [1, 2], "123456", "0", "0", "123456"),
    ]
    model = _make_model_with_rbar(grids, rbar_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    # 6 dependent DOFs at node 2
    assert len(m_set) == 6
    assert GMN.shape == (6, 12)

    # Should be same as RBE2
    m_t1 = m_set[(2, 1)]
    assert_allclose(GMN[m_t1, dof_map[(1, 1)]], 1.0, atol=1e-15)

    m_t2 = m_set[(2, 2)]
    assert_allclose(GMN[m_t2, dof_map[(1, 6)]], 1.0, atol=1e-15)


def test_rbar_global_assembly_split():
    """Test assemble_gmn with RBAR having split independent DOFs."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 2.0, 0.0, 0.0),
    ]
    # CNA=123 (T at A ind), CNB=456 (R at B ind)
    # CMA=456 (R at A dep), CMB=123 (T at B dep)
    rbar_cards = [
        (100, [1, 2], "123", "456", "456", "123"),
    ]
    model = _make_model_with_rbar(grids, rbar_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    assert len(m_set) == 6  # 3 dep at A + 3 dep at B

    # Verify with rigid body
    ndof = 12
    u_full = np.zeros(ndof)
    tx, ty, tz = 0.5, -0.3, 0.1
    rx, ry, rz = 0.01, -0.005, 0.02

    # Set independent DOFs
    u_full[dof_map[(1, 1)]] = tx
    u_full[dof_map[(1, 2)]] = ty
    u_full[dof_map[(1, 3)]] = tz
    u_full[dof_map[(2, 4)]] = rx
    u_full[dof_map[(2, 5)]] = ry
    u_full[dof_map[(2, 6)]] = rz

    u_m = GMN @ u_full

    # Dependent R at A = rx, ry, rz
    assert_allclose(u_m[m_set[(1, 4)]], rx, atol=1e-12)
    assert_allclose(u_m[m_set[(1, 5)]], ry, atol=1e-12)
    assert_allclose(u_m[m_set[(1, 6)]], rz, atol=1e-12)

    # Dependent T at B = translation + omega x offset
    offset_b = np.array([2.0, 0.0, 0.0])
    exp_t_b = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset_b)
    assert_allclose(u_m[m_set[(2, 1)]], exp_t_b[0], atol=1e-12)
    assert_allclose(u_m[m_set[(2, 2)]], exp_t_b[1], atol=1e-12)
    assert_allclose(u_m[m_set[(2, 3)]], exp_t_b[2], atol=1e-12)


# ============================================================================
# Tests: RBAR1 GMN
# ============================================================================


def test_rbar1_coincident():
    """RBAR1 with coincident nodes: GMN is identity for CB DOFs."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 0.0, 0.0, 0.0),
    ]
    rbar1_cards = [(100, [1, 2], "123456")]
    model = _make_model_with_rbar1(grids, rbar1_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    assert len(m_set) == 6
    assert GMN.shape == (6, 12)
    # Identity: each dep DOF at B equals same DOF at A
    for dof in range(1, 7):
        m_idx = m_set[(2, dof)]
        assert_allclose(GMN[m_idx, dof_map[(1, dof)]], 1.0, atol=1e-15)


def test_rbar1_offset():
    """RBAR1 with offset: same as RBE2 with single dependent node."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 2.0, 1.0, 0.0),
    ]
    rbar1_cards = [(100, [1, 2], "123456")]
    model = _make_model_with_rbar1(grids, rbar1_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    # Verify with rigid body motion
    ndof = 12
    u_full = np.zeros(ndof)
    tx, ty, tz = 0.5, -0.3, 0.1
    rx, ry, rz = 0.01, -0.005, 0.02
    for i, val in enumerate([tx, ty, tz, rx, ry, rz]):
        u_full[dof_map[(1, i + 1)]] = val

    u_m = GMN @ u_full

    offset = np.array([2.0, 1.0, 0.0])
    exp_t = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset)
    for i in range(3):
        assert_allclose(u_m[m_set[(2, i + 1)]], exp_t[i], atol=1e-12)
    for i in range(3):
        assert_allclose(u_m[m_set[(2, i + 4)]], [rx, ry, rz][i], atol=1e-12)


def test_rbar1_partial_cb():
    """RBAR1 with CB=123 (only translations dependent)."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
    ]
    rbar1_cards = [(100, [1, 2], "123")]
    model = _make_model_with_rbar1(grids, rbar1_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    assert len(m_set) == 3
    assert GMN.shape == (3, 12)

    # T2 at B: should have coupling from R3 at A (offset x=1)
    m_t2 = m_set[(2, 2)]
    assert_allclose(GMN[m_t2, dof_map[(1, 6)]], 1.0, atol=1e-15)


# ============================================================================
# Tests: RROD GMN
# ============================================================================


def test_rrod_gmn_x_axis_dep_at_b():
    """RROD along x-axis, dependent DOF T1 at node B."""
    xyz_a = np.array([0.0, 0.0, 0.0])
    xyz_b = np.array([1.0, 0.0, 0.0])

    gmn_row, ind_dof_map = compute_rrod_gmn_row(xyz_a, xyz_b, cma=0, cmb=1)

    # Rod along x: e = (1,0,0). dep = uB[0].
    # Constraint: e.(uB - uA) = 0 => uBx = uAx
    # ind DOFs: uA[0], uA[1], uA[2], uB[1], uB[2]
    assert gmn_row.shape == (1, 5)

    # uBx = (e[0]/e[0])*uAx + (e[1]/e[0])*uAy + (e[2]/e[0])*uAz
    #      - (e[1]/e[0])*uBy - (e[2]/e[0])*uBz
    # = 1*uAx + 0 + 0 - 0 - 0 = uAx
    # So coefficient for uAx should be 1.0, all others 0
    assert_allclose(gmn_row[0, 0], 1.0, atol=1e-15)  # uAx
    assert_allclose(gmn_row[0, 1], 0.0, atol=1e-15)  # uAy
    assert_allclose(gmn_row[0, 2], 0.0, atol=1e-15)  # uAz
    assert_allclose(gmn_row[0, 3], 0.0, atol=1e-15)  # uBy
    assert_allclose(gmn_row[0, 4], 0.0, atol=1e-15)  # uBz


def test_rrod_gmn_x_axis_dep_at_a():
    """RROD along x-axis, dependent DOF T1 at node A."""
    xyz_a = np.array([0.0, 0.0, 0.0])
    xyz_b = np.array([1.0, 0.0, 0.0])

    gmn_row, ind_dof_map = compute_rrod_gmn_row(xyz_a, xyz_b, cma=1, cmb=0)

    # uAx = (e[0]/e[0])*uBx + ... = uBx
    # ind DOFs: uA[1], uA[2], uB[0], uB[1], uB[2]
    assert gmn_row.shape == (1, 5)
    assert_allclose(gmn_row[0, 0], 0.0, atol=1e-15)  # uAy
    assert_allclose(gmn_row[0, 1], 0.0, atol=1e-15)  # uAz
    assert_allclose(gmn_row[0, 2], 1.0, atol=1e-15)  # uBx
    assert_allclose(gmn_row[0, 3], 0.0, atol=1e-15)  # uBy
    assert_allclose(gmn_row[0, 4], 0.0, atol=1e-15)  # uBz


def test_rrod_gmn_diagonal():
    """RROD along (1,1,0) direction, dependent T1 at B."""
    xyz_a = np.array([0.0, 0.0, 0.0])
    xyz_b = np.array([1.0, 1.0, 0.0])

    gmn_row, ind_dof_map = compute_rrod_gmn_row(xyz_a, xyz_b, cma=0, cmb=1)

    # e = (1,1,0)/sqrt(2), dep = uB[0], ed = e[0] = 1/sqrt(2)
    # uBx = (e[0]/ed)*uAx + (e[1]/ed)*uAy + (e[2]/ed)*uAz
    #      - (e[1]/ed)*uBy - (e[2]/ed)*uBz
    # = 1*uAx + 1*uAy + 0 - 1*uBy - 0
    assert gmn_row.shape == (1, 5)
    assert_allclose(gmn_row[0, 0], 1.0, atol=1e-12)  # uAx: e[0]/ed = 1
    assert_allclose(gmn_row[0, 1], 1.0, atol=1e-12)  # uAy: e[1]/ed = 1
    assert_allclose(gmn_row[0, 2], 0.0, atol=1e-12)  # uAz: e[2]/ed = 0
    assert_allclose(gmn_row[0, 3], -1.0, atol=1e-12)  # uBy: -e[1]/ed = -1
    assert_allclose(gmn_row[0, 4], 0.0, atol=1e-12)  # uBz: -e[2]/ed = 0


def test_rrod_constraint_satisfied():
    """RROD: verify axial constraint e.(uB-uA)=0 is satisfied."""
    xyz_a = np.array([1.0, 2.0, 3.0])
    xyz_b = np.array([4.0, 5.0, 6.0])
    L_vec = xyz_b - xyz_a
    e = L_vec / np.linalg.norm(L_vec)

    # dep = T2 at B (component 2 = 0-based index 1)
    gmn_row, ind_dof_map = compute_rrod_gmn_row(xyz_a, xyz_b, cma=0, cmb=2)

    # Set random independent DOFs
    rng = np.random.default_rng(42)
    u_ind = rng.standard_normal(gmn_row.shape[1])

    # Compute dependent DOF
    u_dep = float((gmn_row @ u_ind)[0])

    # Reconstruct full displacement
    u_a = np.zeros(3)
    u_b = np.zeros(3)
    for i_col, (node_flag, dof_0) in enumerate(ind_dof_map):
        if node_flag == 0:
            u_a[dof_0] = u_ind[i_col]
        else:
            u_b[dof_0] = u_ind[i_col]
    u_b[1] = u_dep  # T2 at B is the dependent DOF (0-based index 1)

    # Check constraint
    axial_elongation = e @ (u_b - u_a)
    assert_allclose(axial_elongation, 0.0, atol=1e-12)


def test_rrod_global_assembly():
    """Test assemble_gmn with an RROD element."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
    ]
    # Rod along x, dependent T1 at node 2
    rrod_cards = [(100, [1, 2], "", "1")]
    model = _make_model_with_rrod(grids, rrod_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    # 1 dependent DOF: T1 at node 2
    assert len(m_set) == 1
    assert (2, 1) in m_set
    assert GMN.shape == (1, 12)

    # For x-axis rod: uB_x = uA_x (coefficient 1.0)
    m_idx = m_set[(2, 1)]
    assert_allclose(GMN[m_idx, dof_map[(1, 1)]], 1.0, atol=1e-15)
    # No coupling to transverse
    assert_allclose(GMN[m_idx, dof_map[(1, 2)]], 0.0, atol=1e-15)
    assert_allclose(GMN[m_idx, dof_map[(1, 3)]], 0.0, atol=1e-15)


def test_rrod_global_assembly_diagonal():
    """RROD along diagonal: verify constraint via global assembly."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 1.0, 1.0),
    ]
    # dep = T1 at A
    rrod_cards = [(100, [1, 2], "1", "")]
    model = _make_model_with_rrod(grids, rrod_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    assert len(m_set) == 1
    assert (1, 1) in m_set

    # Verify constraint: set some independent DOFs, check e.(uB-uA)=0
    ndof = 12
    u_full = np.zeros(ndof)
    # Set u_A_y=0.5, u_A_z=-0.3, u_B=(0.1, 0.2, 0.4)
    u_full[dof_map[(1, 2)]] = 0.5
    u_full[dof_map[(1, 3)]] = -0.3
    u_full[dof_map[(2, 1)]] = 0.1
    u_full[dof_map[(2, 2)]] = 0.2
    u_full[dof_map[(2, 3)]] = 0.4

    u_m = GMN @ u_full
    u_a1 = u_m[m_set[(1, 1)]]

    # Verify constraint
    e = np.array([1, 1, 1], dtype=float) / np.sqrt(3)
    u_a = np.array([u_a1, 0.5, -0.3])
    u_b = np.array([0.1, 0.2, 0.4])
    assert_allclose(e @ (u_b - u_a), 0.0, atol=1e-12)


# ============================================================================
# Tests: Combined (all types)
# ============================================================================


def test_combined_all_types_assembly():
    """Model with all 6 rigid element types — verify shapes and m-set.

    Layout:
        Nodes 1-2: RBAR (1=ind 123456, 2=dep 123456)
        Nodes 3-4: RBE2 (3=ind 123456, 4=dep)
        Nodes 5-7: RBE3 (5=dep refgrid 123, 6,7=ind)
        Nodes 8-10: RBE1 (8,9=GN ind, 10=GM dep)
        Nodes 11-12: RBAR1 (11=ind, 12=dep CB=123456)
        Nodes 13-14: RROD (along x, dep T1 at 14)
    """
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
        (3, 3.0, 0.0, 0.0),
        (4, 4.0, 0.0, 0.0),
        (5, 6.0, 0.0, 0.0),
        (6, 7.0, 0.0, 0.0),
        (7, 5.0, 0.0, 0.0),
        (8, 10.0, 0.0, 0.0),
        (9, 12.0, 0.0, 0.0),
        (10, 11.0, 1.0, 0.0),
        (11, 15.0, 0.0, 0.0),
        (12, 16.0, 0.0, 0.0),
        (13, 20.0, 0.0, 0.0),
        (14, 21.0, 0.0, 0.0),
    ]
    rbar_cards = [
        (10, [1, 2], "123456", "0", "0", "123456"),
    ]
    rbe2_cards = [
        (20, 3, "123456", [4]),
    ]
    rbe3_cards = [
        (30, 5, "123", [(1.0, "123", [6, 7])]),
    ]
    rbe1_cards = [
        (40, [8, 9], ["123", "456"], [10], ["123456"]),
    ]
    rbar1_cards = [
        (50, [11, 12], "123456"),
    ]
    rrod_cards = [
        (60, [13, 14], "", "1"),
    ]
    model = _make_model_all_types(
        grids,
        rbe1_cards=rbe1_cards,
        rbe2_cards=rbe2_cards,
        rbe3_cards=rbe3_cards,
        rbar_cards=rbar_cards,
        rbar1_cards=rbar1_cards,
        rrod_cards=rrod_cards,
    )

    GMN, m_set, dof_map = assemble_gmn(model)

    # Count m-set DOFs:
    # RBE1: node 10, 6 DOFs
    # RBE2: node 4, 6 DOFs
    # RBE3: node 5, 3 DOFs (comp 123)
    # RBAR: node 2, 6 DOFs
    # RBAR1: node 12, 6 DOFs
    # RROD: node 14, 1 DOF (T1)
    assert len(m_set) == 28
    assert GMN.shape == (28, 84)  # 14 nodes x 6 DOFs = 84


def test_combined_all_types_rigid_body():
    """All 6 element types with consistent rigid-body translation.

    Verify each dependent node receives the same translation when
    independent nodes all move with pure translation.
    """
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
        (3, 3.0, 0.0, 0.0),
        (4, 4.0, 0.0, 0.0),
        (5, 6.0, 0.0, 0.0),
        (6, 7.0, 0.0, 0.0),
        (7, 5.0, 0.0, 0.0),
        (8, 10.0, 0.0, 0.0),
        (9, 12.0, 0.0, 0.0),
        (10, 11.0, 1.0, 0.0),
        (11, 15.0, 0.0, 0.0),
        (12, 16.0, 0.0, 0.0),
        (13, 20.0, 0.0, 0.0),
        (14, 21.0, 0.0, 0.0),
    ]
    rbar_cards = [
        (10, [1, 2], "123456", "0", "0", "123456"),
    ]
    rbe2_cards = [
        (20, 3, "123456", [4]),
    ]
    rbe3_cards = [
        (30, 5, "123", [(1.0, "123", [6, 7])]),
    ]
    rbe1_cards = [
        (40, [8, 9], ["123", "456"], [10], ["123456"]),
    ]
    rbar1_cards = [
        (50, [11, 12], "123456"),
    ]
    rrod_cards = [
        (60, [13, 14], "", "1"),  # RROD along x, dep T1 at 14
    ]
    model = _make_model_all_types(
        grids,
        rbe1_cards=rbe1_cards,
        rbe2_cards=rbe2_cards,
        rbe3_cards=rbe3_cards,
        rbar_cards=rbar_cards,
        rbar1_cards=rbar1_cards,
        rrod_cards=rrod_cards,
    )

    GMN, m_set, dof_map = assemble_gmn(model)

    # Apply pure translation to all nodes (rigid body)
    tx, ty, tz = 1.0, -0.5, 0.3
    ndof = 84
    u_full = np.zeros(ndof)
    for nid in range(1, 15):
        u_full[dof_map[(nid, 1)]] = tx
        u_full[dof_map[(nid, 2)]] = ty
        u_full[dof_map[(nid, 3)]] = tz

    u_m = GMN @ u_full

    # All dependent translations should equal tx, ty, tz
    dep_nodes = [2, 4, 5, 10, 12, 14]
    for nid in dep_nodes:
        for dof, expected in [(1, tx), (2, ty), (3, tz)]:
            if (nid, dof) in m_set:
                assert_allclose(
                    u_m[m_set[(nid, dof)]],
                    expected,
                    atol=1e-12,
                    err_msg=f"node {nid} dof {dof}",
                )


def test_combined_rbe2_rbe3_rigid_body():
    """Combined RBE2+RBE3: rigid body motion preserved through both.

    Layout:
        Node 1 (0,0,0) - RBE2 independent
        Node 2 (2,0,0) - RBE2 dependent, also RBE3 independent
        Node 3 (4,0,0) - RBE3 independent
        Node 4 (3,0,0) - RBE3 dependent (refgrid, between nodes 2 and 3)
    """
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 2.0, 0.0, 0.0),
        (3, 4.0, 0.0, 0.0),
        (4, 3.0, 0.0, 0.0),
    ]
    rbe2_cards = [
        (100, 1, "123456", [2]),
    ]
    rbe3_cards = [
        (200, 4, "123", [(1.0, "123", [2, 3])]),
    ]
    model = _make_model_with_rbe2_and_rbe3(grids, rbe2_cards, rbe3_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    # Apply rigid body translation at node 1
    ndof = 24
    u_full = np.zeros(ndof)
    tx, ty, tz = 1.0, 0.5, -0.2

    # Set node 1 translations
    u_full[dof_map[(1, 1)]] = tx
    u_full[dof_map[(1, 2)]] = ty
    u_full[dof_map[(1, 3)]] = tz

    # Also set nodes 2 and 3 consistently (they are n-set for RBE3)
    u_full[dof_map[(2, 1)]] = tx
    u_full[dof_map[(2, 2)]] = ty
    u_full[dof_map[(2, 3)]] = tz
    u_full[dof_map[(3, 1)]] = tx
    u_full[dof_map[(3, 2)]] = ty
    u_full[dof_map[(3, 3)]] = tz

    u_m = GMN @ u_full

    # Node 2 from RBE2: should match node 1 (pure translation, zero rotation)
    for dof in range(1, 7):
        m_idx = m_set[(2, dof)]
        if dof <= 3:
            assert_allclose(u_m[m_idx], [tx, ty, tz][dof - 1], atol=1e-14)
        else:
            assert_allclose(u_m[m_idx], 0.0, atol=1e-14)

    # Node 4 from RBE3: average of nodes 2 and 3 translations
    for dof in range(1, 4):
        m_idx = m_set[(4, dof)]
        assert_allclose(u_m[m_idx], [tx, ty, tz][dof - 1], atol=1e-12)


def test_combined_rbar_rbe1_rigid_body():
    """RBAR + RBE1 both enforcing the same rigid body.

    Layout:
        Nodes 1-2: RBAR (ind at 1, dep at 2)
        Nodes 3-4-5: RBE1 (GN=3 123456, GM=4,5 T123)
    """
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
        (3, 5.0, 0.0, 0.0),
        (4, 6.0, 1.0, 0.0),
        (5, 4.0, -1.0, 0.0),
    ]
    rbar_cards = [
        (10, [1, 2], "123456", "0", "0", "123456"),
    ]
    rbe1_cards = [
        (20, [3], ["123456"], [4, 5], ["123", "123"]),
    ]
    model = _make_model_all_types(grids, rbe1_cards=rbe1_cards, rbar_cards=rbar_cards)

    GMN, m_set, dof_map = assemble_gmn(model)

    # RBAR: 6 dep DOFs at node 2
    # RBE1: 3+3=6 dep DOFs at nodes 4,5
    assert len(m_set) == 12

    # Apply rigid body with rotation
    tx, ty, tz = 0.5, -0.3, 0.1
    rx, ry, rz = 0.01, -0.005, 0.02
    ndof = 30
    u_full = np.zeros(ndof)

    # Set all independent nodes with rigid body motion
    for nid, xyz in [(1, [0, 0, 0]), (3, [5, 0, 0])]:
        offset = np.array(xyz, dtype=float)
        u_t = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset)
        for i in range(3):
            u_full[dof_map[(nid, i + 1)]] = u_t[i]
        for i in range(3):
            u_full[dof_map[(nid, i + 4)]] = [rx, ry, rz][i]

    u_m = GMN @ u_full

    # Check node 2 (RBAR dep): offset (1,0,0) from node 1
    offset_2 = np.array([1.0, 0.0, 0.0])
    exp_t2 = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset_2)
    for i in range(3):
        assert_allclose(u_m[m_set[(2, i + 1)]], exp_t2[i], atol=1e-12)
    for i in range(3):
        assert_allclose(u_m[m_set[(2, i + 4)]], [rx, ry, rz][i], atol=1e-12)

    # Check nodes 4,5 (RBE1 dep): rigid body from node 3
    for nid, xyz in [(4, [6, 1, 0]), (5, [4, -1, 0])]:
        offset = np.array(xyz, dtype=float)
        exp_t = np.array([tx, ty, tz]) + np.cross([rx, ry, rz], offset)
        for i in range(3):
            assert_allclose(u_m[m_set[(nid, i + 1)]], exp_t[i], atol=1e-12)


def test_combined_no_cross_coupling():
    """Separate rigid elements should not cross-couple."""
    grids = [
        (1, 0.0, 0.0, 0.0),
        (2, 1.0, 0.0, 0.0),
        (3, 10.0, 0.0, 0.0),
        (4, 11.0, 0.0, 0.0),
        (5, 20.0, 0.0, 0.0),
        (6, 21.0, 0.0, 0.0),
        (7, 30.0, 0.0, 0.0),
        (8, 31.0, 0.0, 0.0),
        (9, 32.0, 0.0, 0.0),
    ]
    rbar_cards = [(10, [1, 2], "123456", "0", "0", "123456")]
    rbe2_cards = [(20, 3, "123456", [4])]
    rbe3_cards = [(30, 5, "123", [(1.0, "123", [6])])]
    rbe1_cards = [(40, [7], ["123456"], [8, 9], ["123", "123"])]

    model = _make_model_all_types(
        grids,
        rbe1_cards=rbe1_cards,
        rbe2_cards=rbe2_cards,
        rbe3_cards=rbe3_cards,
        rbar_cards=rbar_cards,
    )

    GMN, m_set, dof_map = assemble_gmn(model)

    # Move only node 1 (RBAR independent)
    ndof = 54
    u_full = np.zeros(ndof)
    u_full[dof_map[(1, 1)]] = 1.0

    u_m = GMN @ u_full

    # Node 2 (RBAR dep) should move
    assert abs(u_m[m_set[(2, 1)]]) > 0.5

    # Nodes from other elements should NOT move
    assert_allclose(u_m[m_set[(4, 1)]], 0.0, atol=1e-15)
    assert_allclose(u_m[m_set[(5, 1)]], 0.0, atol=1e-15)
    assert_allclose(u_m[m_set[(8, 1)]], 0.0, atol=1e-15)


if __name__ == "__main__":
    test_expand_dof_string()
    test_rigid_body_matrix_zero_offset()
    test_rigid_body_matrix_offset()

    test_rbe3_gmn_collinear_2node_translations()
    test_rbe3_gmn_single_node()
    test_rbe3_gmn_single_node_offset()
    test_rbe3_gmn_unequal_weights()
    test_rbe3_gmn_symmetric_4nodes()
    test_rbe3_constraint_matrices()
    test_rbe3_assemble_gmn_simple()
    test_rbe3_thermal_load_zero_alpha()
    test_rbe3_thermal_load_nonzero()
    test_rbe3_assemble_gmn_simple_with_thermal()
    test_rbe3_rigid_body_preservation()
    test_rbe3_global_assembly()

    test_rbe2_gmn_coincident()
    test_rbe2_gmn_offset_x()
    test_rbe2_gmn_offset_xyz()
    test_rbe2_gmn_partial_dofs()
    test_rbe2_gmn_rotation_only()
    test_rbe2_rigid_body_preservation()
    test_rbe2_multiple_dependent_nodes()
    test_rbe2_global_assembly_single()
    test_rbe2_global_assembly_rigid_motion()

    test_rbe1_gmn_single_node_6dof()
    test_rbe1_gmn_two_nodes_split_dofs()
    test_rbe1_gmn_rigid_body_preservation()
    test_rbe1_gmn_partial_dependent()
    test_rbe1_global_assembly()

    test_rbar_gmn_all_ind_at_a()
    test_rbar_gmn_split_dofs()
    test_rbar_gmn_rigid_body_preservation()
    test_rbar_gmn_all_ind_at_b()
    test_rbar_global_assembly()
    test_rbar_global_assembly_split()

    test_rbar1_coincident()
    test_rbar1_offset()
    test_rbar1_partial_cb()

    test_rrod_gmn_x_axis_dep_at_b()
    test_rrod_gmn_x_axis_dep_at_a()
    test_rrod_gmn_diagonal()
    test_rrod_constraint_satisfied()
    test_rrod_global_assembly()
    test_rrod_global_assembly_diagonal()

    test_combined_all_types_assembly()
    test_combined_all_types_rigid_body()
    test_combined_rbe2_rbe3_rigid_body()
    test_combined_rbar_rbe1_rigid_body()
    test_combined_no_cross_coupling()

    print("All GMN tests passed!")
